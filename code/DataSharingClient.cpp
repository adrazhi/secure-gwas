#include "connect.h"
#include "mpc.h"
#include "protocol.h"
#include "util.h"
#include "NTL/ZZ_p.h"
#include "gwasiter.h"

#include <cstdlib>
#include <fstream>
#include <map>
#include <iostream>
#include <sstream>
#include <omp.h>
#include <bits/stdc++.h>
#include <sys/time.h>

using namespace NTL;
using namespace std;

bool send_stream(string data_dir, MPCEnv& mpc, int mode, long start_line, long num_lines) {
  string geno_file = data_dir + "geno.txt";
  string pheno_file = data_dir + "pheno.txt";
  string cov_file = data_dir + "cov.txt";

  bool pheno_flag = (mode == GwasIterator::GP_CODE) || (mode == GwasIterator::GMP_CODE);
  bool missing_flag = (mode == GwasIterator::GM_CODE) || (mode == GwasIterator::GMP_CODE);

  ifstream fin_geno(geno_file.c_str());
  ifstream fin_pheno(pheno_file.c_str());
  ifstream fin_cov(cov_file.c_str());

  if (!fin_geno.is_open()) {
    cout << "Error: could not open " << geno_file << endl;
    return false;
  }

  if (pheno_flag && !fin_pheno.is_open()) {
    cout << "Error: could not open " << pheno_file << endl;
    return false;
  }

  if (!fin_cov.is_open()) {
    cout << "Error: could not open " << cov_file << endl;
    return false;
  }

  long val;
  string line;
  uint32_t lineno = 0;

  // first skip ahead to the desired starting line
  while (lineno < start_line) {
    if (!getline(fin_geno, line)) {
      cout << "Error: data matrix does not have NUM_INDS rows" << endl;
      return false;
    }
    if (pheno_flag) {
      if (!getline(fin_cov, line)) {
        cout << "Error: covariate matrix does not have NUM_INDS rows" << endl;
        return false;
      }
      if (!getline(fin_pheno, line)) {
        cout << "Error: phenotype vector does not have NUM_INDS rows" << endl;
        return false;
      }
    }
    lineno++;
  }

  // now read the desired number of lines
  for (long i = 0; i < num_lines; i++) {
    if (!getline(fin_geno, line)) {
      cout << "Error: data matrix does not have NUM_INDS rows" << endl;
      return false;
    }
    istringstream iss_geno(line);

    if (pheno_flag) {
      Vec<ZZ_p> p;
      p.SetLength(1 + Param::NUM_COVS);

      fin_pheno >> p[0];

      if (!getline(fin_cov, line)) {
        cout << cov_file << " has fewer lines than expected" << endl;
        return false;
      }

      istringstream iss_cov(line);
      for (int j = 0; j < Param::NUM_COVS; j++) {
        iss_cov >> val;
        p[j+1] = ZZ_p(val);
      }

      Vec<ZZ_p> rp;
      mpc.SwitchSeed(1);
      mpc.RandVec(rp, 1 + Param::NUM_COVS);
      mpc.RestoreSeed();

      p -= rp;
      mpc.SendVec(p, 2);
    } else {
      iss_geno >> val;
    }

    Mat<ZZ_p> g;
    Vec<ZZ_p> m;
    if (missing_flag) {
      Init(g, 3, Param::NUM_SNPS);
      Init(m, Param::NUM_SNPS);
    } else {
      Init(g, 1, Param::NUM_SNPS);
    }

    // Read from file
    for (int j = 0; j < Param::NUM_SNPS; j++) {
      string str;
      iss_geno >> str;
      if (str == "NA" || str == "-1") {
        val = -1;
      } else if (str == "0") {
        val = 0;
      } else if (str == "1") {
        val = 1;
      } else if (str == "2") {
        val = 2;
      } else {
        cout << "Error: unknown value in dosage matrix (" << str << ")" << endl;
        return false;
      }

      if (missing_flag) {
        if (val == 0) {
          g[0][j] = 1;
        } else if (val == 1) {
          g[1][j] = 1;
        } else if (val == 2) {
          g[2][j] = 1;
        } else {
          m[j] = 1;
        }
      } else {
        g[0][j] = ZZ_p(val);
      }
    }

    // Generate masks
    Mat<ZZ_p> rg;
    Vec<ZZ_p> rm;
    mpc.SwitchSeed(1);
    if (missing_flag) {
      mpc.RandMat(rg, 3, Param::NUM_SNPS);
      mpc.RandVec(rm, Param::NUM_SNPS);
    } else {
      mpc.RandMat(rg, 1, Param::NUM_SNPS);
    }
    mpc.RestoreSeed();

    // Send masked data
    g -= rg;
    mpc.SendMat(g, 2);

    if (missing_flag) {
      m -= rm;
      mpc.SendVec(m, 2);
    }
  }

  fin_geno.close();
  fin_pheno.close();
  fin_cov.close();

  return true;
}

int main(int argc, char** argv) {
  if (argc < 4) {
    cout << "Usage: DataSharingClient party_id param_file round_number [num_threads] [data_dir (for P3/SP)]" << endl;
    return 1;
  }

  string pid_str(argv[1]);
  int pid;
  if (!Param::Convert(pid_str, pid, "party_id") || pid < 0 || pid > 3) {
    cout << "Error: party_id should be 0, 1, 2, or 3" << endl;
    return 1;
  }
  
  if (!Param::ParseFile(argv[2])) {
    cout << "Could not finish parsing parameter file" << endl;
    return 1;
  }

  string round_str(argv[3]);
  Param::Convert(round_str, Param::CUR_ROUND, "CUR_ROUND");

  // If not in chunked mode, Data Sharing code should not be multi-threaded
  if (!Param::CHUNK_MODE) {
    Param::NUM_THREADS = 1;
  }
  
  string data_dir;
  if (pid == 3) {
    if (argc < 5) {
      cout << "Error: for P3/SP, data directory should be provided as the last argument" << endl;
      return 1;
    }

    data_dir = argv[4];
    if (data_dir[data_dir.size() - 1] != '/') {
      data_dir += "/";
    }

    cout << "Data directory: " << data_dir << endl;
  }

  vector< pair<int, int> > pairs;
  pairs.push_back(make_pair(0, 1));
  pairs.push_back(make_pair(0, 2));
  pairs.push_back(make_pair(1, 2));
  pairs.push_back(make_pair(1, 3));
  pairs.push_back(make_pair(2, 3));

  /* Initialize MPC environment */
  MPCEnv mpc;
  if (!mpc.Initialize(pid, pairs)) {
    cout << "MPC environment initialization failed" << endl;
    return 1;
  }

  // divide up the dataset into chunks
  int num_chunks = 1;
  if (Param::CHUNK_MODE) {
    num_chunks = Param::NUM_CHUNKS[Param::CUR_ROUND];
  }
  long total_lines = Param::NUM_INDS[Param::CUR_ROUND];
  long chunk_size = ceil(total_lines / ((double) num_chunks));
  num_chunks = ceil(total_lines / ((double) chunk_size)); // this is necessary to avoid edge case due to ceiling operator

  // profile runtime for end to end data sharing
  struct timeval start, end;
  double runtime;

  gettimeofday(&start, NULL); 
  ios_base::sync_with_stdio(false);
  bool success = true;
  if (pid < 3) {
    #pragma omp parallel for num_threads(Param::NUM_THREADS)
    for (int i = 0; i < num_chunks; i++) {
      long start_line = chunk_size * i;
      long end_line = start_line + chunk_size;
      if (end_line > total_lines) {
        end_line = total_lines;
      }
      long num_lines = end_line - start_line;

      // avoid error by re-setting the modulus within each thread
      ZZ base_p = conv<ZZ>(Param::BASE_P.c_str());
      ZZ_p::init(base_p);

      success = data_sharing_protocol(mpc, pid, num_lines, i);
    }
  } else {
    #pragma omp parallel for num_threads(Param::NUM_THREADS)
    for (int i = 0; i < num_chunks; i++) {
      long start_line = chunk_size * i;
      long end_line = start_line + chunk_size;
      if (end_line > total_lines) {
        end_line = total_lines;
      }
      long num_lines = end_line - start_line;

      // avoid error by re-setting the modulus within each thread
      ZZ base_p = conv<ZZ>(Param::BASE_P.c_str());
      ZZ_p::init(base_p);

      /* Stream data upon request */
      int signal = mpc.ReceiveInt(1);

      while (signal != GwasIterator::TERM_CODE) {
        success = send_stream(data_dir, mpc, signal, start_line, num_lines);
        signal = mpc.ReceiveInt(1);
      }
    }

    cout << "Done with streaming data" << endl;
  }
  gettimeofday(&end, NULL);

  // calculate runtime in seconds
  runtime = (end.tv_sec - start.tv_sec) * 1e6;
  runtime = (runtime + (end.tv_usec - start.tv_usec)) * 1e-6;

  // This is here just to keep P0 online until the end for data transfer
  // In practice, P0 would send data in advance before each phase and go offline
  if (pid == 0) {
    mpc.ReceiveBool(2);
  } else if (pid == 2) {
    mpc.SendBool(true, 0);
  }

  mpc.CleanUp();

  if (success) {
    cout << "Protocol successfully completed" << endl;
    cout << "Data Sharing Runtime: " << fixed << runtime << setprecision(6); 
    cout << " sec" << endl;
    return 0;
  } else {
    cout << "Protocol abnormally terminated" << endl;
    return 1;
  }
}
