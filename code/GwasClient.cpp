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
#include <bits/stdc++.h>
#include <sys/time.h>

using namespace NTL;
using namespace std;

void print_vec(const char* name, vector<double> &vect, int num) {
  cout << name << ": ";
  for (int i = 0; i < num; i++) {
    cout << vect[i] << "  ";
  }
  cout << endl;
}

int main(int argc, char** argv) {
  if (argc < 3) {
    cout << "Usage: GwasClient party_id param_file" << endl;
    return 1;
  }

  string pid_str(argv[1]);
  int pid;
  if (!Param::Convert(pid_str, pid, "party_id") || pid < 0 || pid > 2) {
    cout << "Error: party_id should be 0, 1, or 2" << endl;
    return 1;
  }

  // reset the current round
  Param::CUR_ROUND = 0;

  if (!Param::ParseFile(argv[2])) {
    cout << "Could not finish parsing parameter file" << endl;
    return 1;
  }

  // pre-process the param before running GWAS
  int n = Param::NUM_INDS.size();
  for (int i = 0; i < n; i++) {
    int num_chunks = Param::NUM_CHUNKS[i];

    // expand out the numbers of individuals
    long total = Param::NUM_INDS[i];
    long chunk_size = ceil(total / ((double) num_chunks));
    num_chunks = ceil(total / ((double) chunk_size)); 
    long remainder = total - ((num_chunks - 1) * chunk_size);
    for (int j = 0; j < num_chunks - 1; j++) {
      Param::NUM_INDS.push_back(chunk_size);
    }
    Param::NUM_INDS.push_back(remainder);

    // expand out the cache file prefixes
    for (int j = 0; j < num_chunks; j++) {
      Param::CACHE_FILE_PREFIX.push_back(Param::CACHE_FILE_PREFIX[i] + "_" + to_string(j));
    }
  }
  for (int i = 0; i < n; i++) {
    Param::NUM_INDS.pop_front();
    Param::CACHE_FILE_PREFIX.pop_front();
  }
  print_vec("NUM_INDS", Param::NUM_INDS, Param::NUM_INDS.size());
  print_vec("CACHE_FILE_PREFIX", Param::CACHE_FILE_PREFIX, Param::CACHE_FILE_PREFIX.size());

  if (argc == 4) {
    string num_threads_str(argv[3]);
    int num_threads = stoi(num_threads_str);
    Param::NUM_THREADS = num_threads;
  }
  Param::NUM_THREADS = 1; // comment this out once we make GWAS parallel
  cout << "Number of threads: " << Param::NUM_THREADS << endl;

  vector< pair<int, int> > pairs;
  pairs.push_back(make_pair(0, 1));
  pairs.push_back(make_pair(0, 2));
  pairs.push_back(make_pair(1, 2));

  /* Initialize MPC environment */
  MPCEnv mpc;
  if (!mpc.Initialize(pid, pairs)) {
    cout << "MPC environment initialization failed" << endl;
    return 1;
  }

  struct timeval start, end;
  double runtime;

  gettimeofday(&start, NULL); 
  ios_base::sync_with_stdio(false);
  bool success = gwas_protocol(mpc, pid);
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
    cout << "GWAS Runtime: " << fixed << runtime << setprecision(6); 
    cout << " sec" << endl;
    return 0;
  } else {
    cout << "Protocol abnormally terminated" << endl;
    return 1;
  }
}
