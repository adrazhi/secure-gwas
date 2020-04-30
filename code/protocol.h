#ifndef __PROTOCOL_H_
#define __PROTOCOL_H_

#include "gwasiter.h"
#include "mpc.h"
#include "util.h"
#include <vector>
#include <NTL/mat_ZZ_p.h>
#include <NTL/mat_ZZ.h>
#include <NTL/ZZ.h>
#include <NTL/BasicThreadPool.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

#include <chrono>

using namespace NTL;
using namespace std;

using msec = chrono::milliseconds;
using get_time = chrono::steady_clock;

#define ABS(a) (((a)<0)?-(a):(a))

auto clock_start = get_time::now();

void tic() {
  clock_start = get_time::now();
}

int toc() {
  auto clock_end = get_time::now();
  int duration = chrono::duration_cast<msec>(clock_end - clock_start).count();
  cout << "Elapsed time is " << duration / 1000.0 << " secs" << endl;
  return duration;
}

bool DecrComp(const pair<int, double> &a, const pair<int, double> &b) {
  return a.second > b.second;
}

string cache(int pid, string desc) {
  ostringstream oss;
  oss << Param::CACHE_FILE_PREFIX[Param::CUR_ROUND] << "_" << desc << ".bin";
  return oss.str();
}

string cache(int pid, int index, string desc) {
  ostringstream oss;
  oss << Param::CACHE_FILE_PREFIX[index] << "_" << desc << ".bin";
  return oss.str();
}

string cache(int pid, int index, int chunk_id, string desc) {
  ostringstream oss;
  oss << Param::CACHE_FILE_PREFIX[index] << "_" << chunk_id << "_" << desc << ".bin";
  return oss.str();
}

string cache(int pid, int index) {
  ostringstream oss;
  oss << Param::CACHE_FILE_PREFIX[Param::CUR_ROUND] << "_" << index << ".bin";
  return oss.str();
}

string outname(string desc) {
  ostringstream oss;
  oss << Param::OUTPUT_FILE_PREFIX << "_" << desc << ".txt";
  return oss.str();
}

bool logireg_protocol(MPCEnv& mpc, int pid) {
  SetNumThreads(Param::NUM_THREADS);
  cout << AvailableThreads() << " threads created" << endl;

  int ntop = 100;

  // need to update logireg_protocol function to be compatible with multiple input datasets
  int n0 = Param::NUM_INDS[Param::CUR_ROUND];
  int m0 = Param::NUM_SNPS;
  int k = Param::NUM_DIM_TO_REMOVE;

  cout << "n0: " << n0 << ", " << "m0: " << m0 << endl;

  // Shared variables
  string s;
  fstream fs;
  ofstream ofs;
  ifstream ifs;
  streampos strpos;
  int ind;
  Vec<ZZ_p> tmp_vec;
  Mat<ZZ_p> tmp_mat;

  //mpc.ProfilerPushState("main");

  ZZ_p fp_one = DoubleToFP(1, Param::NBIT_K, Param::NBIT_F);

  Vec<ZZ_p> pheno;
  Init(pheno, n0);

  Mat<ZZ_p> cov;
  Init(cov, n0, Param::NUM_COVS);

  if (!exists(cache(pid, "input_geno")) || !exists(cache(pid, "input_pheno_cov"))) {
    cout << "Initial data sharing results not found:" << endl;
    cout << "\t" << cache(pid, "input_geno") << endl;
    cout << "\t" << cache(pid, "input_pheno_cov") << endl;
    return false;
  }

  cout << "Initial data sharing results found" << endl;

  ifs.open(cache(pid, "input_pheno_cov").c_str(), ios::binary);
  mpc.ReadFromFile(pheno, ifs, n0);
  mpc.ReadFromFile(cov, ifs, n0, Param::NUM_COVS);
  ifs.close();

  cout << "Phenotypes and covariates loaded" << endl;

  Vec<ZZ_p> gkeep1;
  Init(gkeep1, m0);
  
  cout << "Using locus missing rate filter from a previous run" << endl;
  
  if (pid == 2) {
    ifs.open(outname("gkeep1").c_str());
    for (int i = 0; i < m0; i++) {
      ifs >> gkeep1[i];
    }
    ifs.close();

    mpc.SendVec(gkeep1, 0);
    mpc.SendVec(gkeep1, 1);
  } else {
    mpc.ReceiveVec(gkeep1, 2, m0);
  }

  uint m1 = conv<uint>(Sum(gkeep1));
  cout << "n0: " << n0 << ", " << "m1: " << m1 << endl;

  Vec<ZZ_p> ikeep;
  Init(ikeep, n0);

  cout << "Using individual missing rate/het rate filters from a previous run" << endl;
  
  if (pid == 2) {
    ifs.open(outname("ikeep").c_str());
    for (int i = 0; i < n0; i++) {
      ifs >> ikeep[i];
    }
    ifs.close();

    mpc.SendVec(ikeep, 0);
    mpc.SendVec(ikeep, 1);
  } else {
    mpc.ReceiveVec(ikeep, 2, n0);
  }

  uint n1 = conv<uint>(Sum(ikeep));
  cout << "n1: " << n1 << ", " << "m1: " << m1 << endl;

  cout << "Filtering phenotypes and covariates" << endl;
  mpc.Filter(pheno, ikeep, n1);
  mpc.FilterRows(cov, ikeep, n1);

  Vec<ZZ_p> gkeep2;
  Init(gkeep2, m1);

  cout << "Using MAF/HWE filters from a previous run" << endl;
  
  if (pid == 2) {
    ifs.open(outname("gkeep2").c_str());
    for (int i = 0; i < m1; i++) {
      ifs >> gkeep2[i];
    }
    ifs.close();

    mpc.SendVec(gkeep2, 0);
    mpc.SendVec(gkeep2, 1);
  } else {
    mpc.ReceiveVec(gkeep2, 2, m1);
  }

  uint m2 = conv<uint>(Sum(gkeep2));
  cout << "n1: " << n1 << ", " << "m2: " << m2 << endl;

  cout << "Using CA statistics from a previous run" << endl;
  
  Vec<ZZ_p> gkeep3;
  Init(gkeep3, m2);

  if (pid == 2) {
    vector<pair<int, double> > cavec(m2);

    ifs.open(outname("assoc").c_str());
    double val;
    for (int i = 0; i < m2; i++) {
      ifs >> val;
      cavec[i] = make_pair(i, val * val);
    }
    ifs.close();

    sort(cavec.begin(), cavec.end(), DecrComp);

    cout << "Selected top " << ntop << " candidates" << endl;
    cout << "Top 5 CA stats: " << cavec[0].second;
    for (int i = 1; i < 5; i++) {
      cout << ", " << cavec[i].second;
    }
    cout << endl;

    for (int i = 0; i < ntop; i++) {
      gkeep3[cavec[i].first] = 1;
    }
    mpc.SendVec(gkeep3, 0);
    mpc.SendVec(gkeep3, 1);
  } else {
    mpc.ReceiveVec(gkeep3, 2, m2);
  }

  Mat<ZZ_p> V;
  Init(V, k, n1);

  cout << "Using eigenvectors from a previous run" << endl;
  ifs.open(cache(pid, "eigen").c_str(), ios::binary);
  mpc.ReadFromFile(V, ifs, k, n1);
  ifs.close();

  // Concatenate covariate matrix and jointly orthogonalize
  mpc.Transpose(cov);
  V.SetDims(k + Param::NUM_COVS, n1);
  if (pid > 0) {
    for (int i = 0; i < Param::NUM_COVS; i++) {
      V[k + i] = cov[i] * fp_one;
    }
  }
  cov.kill();

  mpc.OrthonormalBasis(V, V);

  Vec<ZZ_p> V_mean;
  Init(V_mean, V.NumRows());
  ZZ_p fp_denom = DoubleToFP(1 / ((double) V.NumCols()), Param::NBIT_K, Param::NBIT_F);
  for (int i = 0; i < V_mean.length(); i++) {
    V_mean[i] = Sum(V[i]) * fp_denom;
  }
  mpc.Trunc(V_mean);
  for (int i = 0; i < V_mean.length(); i++) {
    AddScalar(V[i], -V_mean[i]);
  }

  Vec<ZZ_p> V_var;
  mpc.InnerProd(V_var, V);
  mpc.Trunc(V_var);
  V_var *= fp_denom;
  mpc.Trunc(V_var);

  Vec<ZZ_p> V_stdinv, dummy_vec;
  mpc.FPSqrt(dummy_vec, V_stdinv, V_var);
  
  for (int i = 0; i < V_mean.length(); i++) {
    mpc.MultMat(V[i], V[i], V_stdinv[i]);
  }
  mpc.Trunc(V);

  Vec<bool> gkeep;
  gkeep.SetLength(m0);
  for (int j = 0; j < m0; j++) {
    gkeep[j] = gkeep1[j] == 1;
  }

  ind = 0;
  for (int j = 0; j < m0; j++) {
    if (gkeep[j]) {
      gkeep[j] = gkeep2[ind] == 1;
      ind++;
    }
  }

  ind = 0;
  for (int j = 0; j < m0; j++) {
    if (gkeep[j]) {
      gkeep[j] = gkeep3[ind] == 1;
      ind++;
    }
  }

  Mat<ZZ_p> X, X_mask;
  if (exists(cache(pid, "logi_input"))) {
    cout << "logi_input cache found" << endl;
    ifs.open(cache(pid, "logi_input").c_str(), ios::in | ios::binary);
    mpc.BeaverReadFromFile(X, X_mask, ifs, ntop, n1);
    ifs.close();
  } else {
    X.SetDims(n1, ntop);
    X_mask.SetDims(n1, ntop);

    ifs.open(cache(pid, "input_geno").c_str(), ios::binary);
    if (pid > 0) {
      mpc.ImportSeed(10, ifs);
    } else {
      for (int p = 1; p <= 2; p++) {
        mpc.ImportSeed(10 + p, ifs);
      }
    }

    cout << "Collecting genotypes for top candidates" << endl;
    ind = -1;
    int batch_size = n1 / 10;
    tic();
    for (int cur = 0; cur < n1; cur++) {
      ind++;

      if ((cur + 1) % batch_size == 0 || cur == n1 - 1) {
        cout << cur + 1 << "/" << n1 << ", "; toc();
        tic();
      }

      Mat<ZZ_p> g0, g0_mask;
      Vec<ZZ_p> miss0, miss0_mask;

      while (ikeep[ind] != 1) {
        if (pid > 0) {
          mpc.SkipData(ifs, 3, m0); // g
          mpc.SkipData(ifs, m0); // miss

          mpc.SwitchSeed(10);
          mpc.RandMat(g0_mask, 3, m0);
          mpc.RandVec(miss0_mask, m0);
          mpc.RestoreSeed();
        } else {
          for (int p = 1; p <= 2; p++) {
            mpc.SwitchSeed(10 + p);
            mpc.RandMat(g0_mask, 3, m0);
            mpc.RandVec(miss0_mask, m0);
            mpc.RestoreSeed();
          }
        }
        ind++;
      }

      if (pid > 0) {
        mpc.ReadFromFile(g0, ifs, 3, m0); // g
        mpc.ReadFromFile(miss0, ifs, m0); // miss

        mpc.SwitchSeed(10);
        mpc.RandMat(g0_mask, 3, m0);
        mpc.RandVec(miss0_mask, m0);
        mpc.RestoreSeed();
      } else {
        Init(g0, 3, m0);
        Init(g0_mask, 3, m0);
        Init(miss0, m0);
        Init(miss0_mask, m0);

        for (int p = 1; p <= 2; p++) {
          mpc.SwitchSeed(10 + p);
          mpc.RandMat(tmp_mat, 3, m0);
          mpc.RandVec(tmp_vec, m0);
          mpc.RestoreSeed();

          g0_mask += tmp_mat;
          miss0_mask += tmp_vec;
        }
      }
      
      Mat<ZZ_p> g, g_mask;
      Vec<ZZ_p> miss, miss_mask;
      g.SetDims(3, ntop);
      miss.SetLength(ntop);
      g_mask.SetDims(3, ntop);
      miss_mask.SetLength(ntop);
      int ind2 = 0;
      for (int j = 0; j < m0; j++) {
        if (gkeep[j]) {
          for (int k = 0; k < 3; k++) {
            g[k][ind2] = g0[k][j];
            g_mask[k][ind2] = g0_mask[k][j];
          }
          miss[ind2] = miss0[j];
          miss_mask[ind2] = miss0_mask[j];
          ind2++;
        }
      }

      X[cur] = g[1] + 2 * g[2];
      X_mask[cur] = g_mask[1] + 2 * g_mask[2];
    }

    mpc.Transpose(X); // ntop-by-n1
    transpose(X_mask, X_mask);

    fs.open(cache(pid, "logi_input").c_str(), ios::out | ios::binary);
    mpc.BeaverWriteToFile(X, X_mask, fs);
    fs.close();
  }

  // Shuffle
  Mat<ZZ_p> xrt, xmt;
  transpose(xrt, X);
  transpose(xmt, X_mask);

  transpose(V, V);

  Vec<long> indices;
  indices.SetLength(n1);
  for (int i = 0; i < n1; i++) {
    indices[i] = i + 1;
  }

  mpc.SwitchSeed(-1);
  for (int i = 0; i < n1-1; i++) {
    long chosen = mpc.RandBnd(n1-i);
    Vec<ZZ_p> tmp;
    if (chosen > 0) {
      tmp = xrt[i];
      xrt[i] = xrt[i + chosen];
      xrt[i + chosen] = tmp;

      tmp = xmt[i];
      xmt[i] = xmt[i + chosen];
      xmt[i + chosen] = tmp;

      tmp = V[i];
      V[i] = V[i + chosen];
      V[i + chosen] = tmp;

      ZZ_p tmpy;
      tmpy = pheno[i];
      pheno[i] = pheno[i + chosen];
      pheno[i + chosen] = tmpy;

      long t = indices[i];
      indices[i] = indices[i + chosen];
      indices[i + chosen] = t;
    }
  }
  mpc.RestoreSeed();

  transpose(X, xrt);
  transpose(X_mask, xmt);
  transpose(V, V);
  xrt.kill();
  xmt.kill();

  Mat<ZZ_p> V_mask;
  mpc.BeaverPartition(V_mask, V);

  Vec<ZZ_p> pheno_mask;
  mpc.BeaverPartition(pheno_mask, pheno);

  Vec<ZZ_p> b0;
  Mat<ZZ_p> bv;
  Vec<ZZ_p> bx;
  mpc.ParallelLogisticRegression(b0, bv, bx, X, X_mask, V, V_mask, pheno, pheno_mask, 500);

  fs.open(cache(pid, "logireg_final_coeff").c_str(), ios::out | ios::binary);
  if (pid > 0) {
    mpc.WriteToFile(b0, fs);
    mpc.WriteToFile(bv, fs);
    mpc.WriteToFile(bx, fs);
  }
  fs.close();

  mpc.RevealSym(bx);
  if (pid == 2) {
    Vec<double> bx_double;
    FPToDouble(bx_double, bx, Param::NBIT_K, Param::NBIT_F);
    ofs.open(outname("logi_coeff").c_str(), ios::out);
    for (int i = 0; i < bx_double.length(); i++) {
      ofs << bx_double[i] << endl;
    }
    ofs.close();
    cout << "Result written to " << outname("logi_coeff") << endl;
  }

  return true;
}

bool data_sharing_protocol(MPCEnv& mpc, int pid, int n, int chunk_id) {
  cout << "n: " << n << endl;

  fstream fs;

  Vec<ZZ_p> pheno;
  Init(pheno, n);

  Mat<ZZ_p> cov;
  Init(cov, n, Param::NUM_COVS);

  fs.open(cache(pid, Param::CUR_ROUND, chunk_id, "input_geno").c_str(), ios::out | ios::binary);
  if (pid > 0) {
    mpc.ExportSeed(fs, 0);
  } else {
    for (int p = 1; p <= 2; p++) {
      mpc.ExportSeed(fs, p);
    }
  }

  GwasIterator git(mpc, pid);

  git.Init(true, true);

  long bsize = n / 10;

  cout << "Begin processing:" << endl;

  tic();
  for (int i = 0; i < n; i++) {
    Mat<ZZ_p> g;
    Vec<ZZ_p> miss, p;

    git.GetNextGMP(g, miss, p);

    if (pid > 0) {
      pheno[i] = p[0];
      for (int j = 0; j < Param::NUM_COVS; j++) {
        cov[i][j] = p[1 + j];
      }
    }

    // In practice this would be done in one batch
    Mat<ZZ_p> g_mask;
    Vec<ZZ_p> miss_mask;
    mpc.BeaverPartition(g_mask, g);
    mpc.BeaverPartition(miss_mask, miss);

    if (pid > 0) {
      // Note: g_mask and miss_mask can be recovered from PRG and
      // need not be written
      mpc.WriteToFile(g, fs);
      mpc.WriteToFile(miss, fs);
    }

    if ((i + 1) % bsize == 0 || i == n - 1) {
      cout << "\t" << i+1 << " / " << n << ", "; toc(); tic();
    }
  }

  git.Terminate();

  fs.close();

  cout << "Finished writing Beaver partitioned genotype data" << endl;

  if (Param::DEBUG) {
    cout << "pheno" << endl;
    mpc.Print(pheno, 5);
    cout << "cov" << endl;
    mpc.Print(cov[0], 5);
  }

  fs.open(cache(pid, Param::CUR_ROUND, chunk_id, "input_pheno_cov").c_str(), ios::out | ios::binary);
  mpc.WriteToFile(pheno, fs);
  mpc.WriteToFile(cov, fs);
  fs.close();

  cout << "Finished writing phenotype and covariate data" << endl;

  return true;
}

bool gwas_protocol(MPCEnv& mpc, int pid) {
  SetNumThreads(Param::NTL_NUM_THREADS);
  cout << AvailableThreads() << " threads created for NTL" << endl;
  
  int n0 = 0; // total number of individuals across datasets
  for (int i = 0; i < Param::NUM_INDS.size(); i++) {
    n0 += Param::NUM_INDS[i];
  }
  int m0 = Param::NUM_SNPS;
  int k = Param::NUM_DIM_TO_REMOVE;
  int kp = k + Param::NUM_OVERSAMPLE;

  cout << "n0: " << n0 << ", " << "m0: " << m0 << endl;

  // Shared variables
  string s;
  fstream fs;
  ofstream ofs;
  ifstream ifs;
  streampos strpos;
  int ind;
  Vec<ZZ_p> tmp_vec;
  Mat<ZZ_p> tmp_mat;

  int num_datasets = Param::NUM_INDS.size();
  int num_threads = Param::NUM_THREADS;

  mpc.ProfilerPushState("main");

  // Read in SNP list
  Vec<ZZ> snp_pos;
  Init(snp_pos, m0);

  ifs.open(Param::SNP_POS_FILE.c_str());
  if (!ifs.is_open()) {
    cout << "Could not open SNP_POS_FILE: " << Param::SNP_POS_FILE << endl;
    return false;
  }

  for (int i = 0; i < m0; i++) {
    long chrom, pos;
    ifs >> chrom >> pos;
    snp_pos[i] = ZZ(chrom) * 1e9 + ZZ(pos);
  }
  ifs.close();

  /* Useful constants */
  ZZ_p two(2), twoinv;
  inv(twoinv, two);

  ZZ_p fp_one = DoubleToFP(1, Param::NBIT_K, Param::NBIT_F);

  Vec<ZZ_p> pheno;
  Init(pheno, n0);

  Mat<ZZ_p> cov;
  Init(cov, n0, Param::NUM_COVS);

  for (int i = 0; i < num_datasets; i++) {
    if (!exists(cache(pid, i, "input_geno")) || !exists(cache(pid, i, "input_pheno_cov"))) {
      cout << "Initial data sharing results not found:" << endl;
      cout << "\t" << cache(pid, i, "input_geno") << endl;
      cout << "\t" << cache(pid, i, "input_pheno_cov") << endl;
      return false;
    }
  }

  cout << "Initial data sharing results found" << endl;

  #pragma omp parallel for num_threads(num_threads)
  for (int i = 0; i < num_datasets; i++) {
    long offset = 0;
    for (int j = 0; j < i; j++) {
      offset += Param::NUM_INDS[j];
    }
    long inner_n0 = Param::NUM_INDS[i];

    // avoid error by re-setting the modulus within each thread
    ZZ base_p = conv<ZZ>(Param::BASE_P.c_str());
    ZZ_p::init(base_p);

    ifstream inner_ifs;
    inner_ifs.open(cache(pid, i, "input_pheno_cov").c_str(), ios::binary);

    Vec<ZZ_p> sub_pheno;
    Init(sub_pheno, inner_n0);
    mpc.ReadFromFile(sub_pheno, inner_ifs, inner_n0);
    for (int j = 0; j < inner_n0; j++) {
      pheno[offset + j] = sub_pheno[j];
    }

    Mat<ZZ_p> sub_cov;
    Init(sub_cov, inner_n0, Param::NUM_COVS);
    mpc.ReadFromFile(sub_cov, inner_ifs, inner_n0, Param::NUM_COVS);
    for (int j = 0; j < inner_n0; j++) {
      cov[offset + j] = sub_cov[j];
    }

    inner_ifs.close();
  }

  cout << "Phenotypes and covariates loaded" << endl;

  if (Param::DEBUG) {
    cout << "pheno" << endl;
    mpc.Print(pheno, 5);
    cout << "cov" << endl;
    mpc.Print(cov[0], 5);
  }

  mpc.ProfilerPushState("qc");

  mpc.ProfilerPushState("snp_miss");

  Vec<ZZ_p> gkeep1;
  Init(gkeep1, m0);

  if (Param::SKIP_QC) {
    for (int i = 0; i < m0; i++) {
      gkeep1[i] = 1;
    }
    cout << "Locus missing rate filter skipped" << endl;
  } else {

    bool history;
    if (pid == 2) {
      history = exists(outname("gkeep1"));
      mpc.SendBool(history, 0);
      mpc.SendBool(history, 1);
    } else {
      // ask P2 if gkeep1 has been computed before
      history = mpc.ReceiveBool(2);
    }

    if (history) {
      cout << "Using locus missing rate filter from a previous run" << endl;
   
      if (pid == 2) {
        ifs.open(outname("gkeep1").c_str());
        for (int i = 0; i < m0; i++) {
          ifs >> gkeep1[i];
        }
        ifs.close();

        mpc.SendVec(gkeep1, 0);
        mpc.SendVec(gkeep1, 1);
      } else {
        mpc.ReceiveVec(gkeep1, 2, m0);
      }
    } else {
      Vec<ZZ_p> gmiss;
      Init(gmiss, m0);
      
      if (exists(cache(pid, "gmiss"))) {

        cout << "Locus missing rate cache found" << endl;

        ifs.open(cache(pid, "gmiss").c_str(), ios::binary);
        mpc.ReadFromFile(gmiss, ifs, m0);
        ifs.close();

      } else {

        cout << "Taking a pass to calculate locus missing rates:" << endl;

        if (pid > 0) {
          // Loop over all datasets to calculate missing rate across all individuals
          #pragma omp parallel for num_threads(num_threads)
          for (int i = 0; i < num_datasets; i++) {
            long inner_n0 = Param::NUM_INDS[i];

            // keep track of intermediate results to minimize number of (costly) atomic updates
            Vec<ZZ_p> inner_gmiss;
            Init(inner_gmiss, m0);

            // avoid error by re-setting the modulus within each thread
            ZZ base_p = conv<ZZ>(Param::BASE_P.c_str());
            ZZ_p::init(base_p);

            ifstream inner_ifs;
            inner_ifs.open(cache(pid, i, "input_geno").c_str(), ios::binary);

            mpc.ImportSeed(10, inner_ifs);

            long bsize = inner_n0 / 10;

            tic();
            for (int j = 0; j < inner_n0; j++) {
              Vec<ZZ_p> miss, miss_mask;
              Mat<ZZ_p> skip_mat;

              // Load stored Beaver partition
              mpc.SwitchSeed(10);
              mpc.RandMat(skip_mat, 3, m0); // g_mask
              mpc.RandVec(miss_mask, m0);
              mpc.RestoreSeed();

              if (pid == 2) {
                mpc.SkipData(inner_ifs, 3, m0);
                mpc.ReadFromFile(miss, inner_ifs, m0);
              }

              // Recover secret shares from Beaver partition
              if (pid == 1) {
                miss = miss_mask;
              } else {
                miss += miss_mask;
              }

              inner_gmiss += miss;

              if ((j + 1) % bsize == 0 || j == inner_n0 - 1) {
                cout << "\t" << j+1 << " / " << inner_n0 << ", "; toc(); tic();
              }
            }
            inner_ifs.close();

            // Update global gmiss - this operation must be atomic
            #pragma omp critical ( gmiss_update )
              gmiss += inner_gmiss;
          }
        }

        fs.open(cache(pid, "gmiss").c_str(), ios::out | ios::binary);
        mpc.WriteToFile(gmiss, fs);
        fs.close();

        cout << "Wrote results to cache" << endl;

      }

      if (Param::DEBUG) {
        cout << "gmiss" << endl;
        mpc.Print(gmiss, 5);
      }
      cout << "Locus missing rate filter ... " << endl; tic();

      ZZ_p gmiss_ub = ZZ_p((long) (n0 * Param::GMISS_UB));
      mpc.LessThanPublic(gkeep1, gmiss, gmiss_ub);
      cout << "done. "; toc();
      
      mpc.RevealSym(gkeep1);

      if (pid == 2) {
        mpc.SendVec(gkeep1, 0);
      } else if (pid == 0) {
        mpc.ReceiveVec(gkeep1, 2, m0);
      }

      if (pid == 2) {
        ofs.open(outname("gkeep1").c_str());
        for (int i = 0; i < gkeep1.length(); i++) {
          ofs << gkeep1[i] << endl;
        }
        ofs.close();
      }

    }
  }

  uint m1 = conv<uint>(Sum(gkeep1));
  cout << "n0: " << n0 << ", " << "m1: " << m1 << endl;

  cout << "Filtering SNP position vector ... " << endl; tic();
  FilterVec(snp_pos, gkeep1);
  cout << "done. "; toc();

  mpc.ProfilerPopState(true); // snp_miss

  mpc.ProfilerPushState("ind_miss/het");

  Vec<ZZ_p> ikeep;
  Init(ikeep, n0);

  if (Param::SKIP_QC) {
    for (int i = 0; i < n0; i++) {
      ikeep[i] = 1;
    }
    cout << "Individual missing rate/het rate filters skipped" << endl;
  } else {

    bool history;
    if (pid == 2) {
      history = exists(outname("ikeep"));
      mpc.SendBool(history, 0);
      mpc.SendBool(history, 1);
    } else {
      // ask P2 if ikeep has been computed before
      history = mpc.ReceiveBool(2);
    }

    if (history) {
      cout << "Using individual missing rate/het rate filters from a previous run" << endl;
   
      if (pid == 2) {
        ifs.open(outname("ikeep").c_str());
        for (int i = 0; i < n0; i++) {
          ifs >> ikeep[i];
        }
        ifs.close();

        mpc.SendVec(ikeep, 0);
        mpc.SendVec(ikeep, 1);
      } else {
        mpc.ReceiveVec(ikeep, 2, n0);
      }
    } else {

      Vec<ZZ_p> imiss, ihet;
      Init(imiss, n0);
      Init(ihet, n0);

      if (exists(cache(pid, "imiss_ihet"))) {
        cout << "Individual missing rate and het rate cache found" << endl;

        ifs.open(cache(pid, "imiss_ihet").c_str(), ios::binary);
        mpc.ReadFromFile(imiss, ifs, n0);
        mpc.ReadFromFile(ihet, ifs, n0);
        ifs.close();

      } else {

        cout << "Taking a pass to calculate individual missing rates and het rates:" << endl;

        mpc.ProfilerPushState("data_scan");

        if (pid > 0) {
          // Loop over all datasets to calculate individual missing/heterozygous rate for all individuals
          #pragma omp parallel for num_threads(num_threads)
          for (int i = 0; i < num_datasets; i++) {
            long offset = 0;
            for (int j = 0; j < i; j++) {
              offset += Param::NUM_INDS[j];
            }
            long inner_n0 = Param::NUM_INDS[i];

            // avoid error by re-setting the modulus within each thread
            ZZ base_p = conv<ZZ>(Param::BASE_P.c_str());
            ZZ_p::init(base_p);

            ifstream inner_ifs;
            inner_ifs.open(cache(pid, i, "input_geno").c_str(), ios::binary);

            mpc.ImportSeed(10, inner_ifs);

            long bsize = inner_n0 / 10;

            tic();
            for (int j = 0; j < inner_n0; j++) {
              Mat<ZZ_p> g, g_mask;
              Vec<ZZ_p> miss, miss_mask;

              // Load stored Beaver partition
              mpc.SwitchSeed(10);
              mpc.RandMat(g_mask, 3, m0);
              mpc.RandVec(miss_mask, m0);
              mpc.RestoreSeed();

              if (pid == 2) {
                mpc.ReadFromFile(g, inner_ifs, 3, m0);
                mpc.ReadFromFile(miss, inner_ifs, m0);
              }

              // Recover secret shares from Beaver partition
              if (pid == 1) {
                g = g_mask;
                miss = miss_mask;
              } else {
                g += g_mask;
                miss += miss_mask;
              }

              // Add to running sum
              for (int k = 0; k < m0; k++) {
                if (gkeep1[k] == 1) {
                  imiss[offset + j] += miss[k];
                  ihet[offset + j] += g[1][k];
                }
              }

              if ((j + 1) % bsize == 0 || j == inner_n0 - 1) {
                cout << "\t" << j+1 << " / " << inner_n0 << ", "; toc(); tic();
              }
            }
            inner_ifs.close();
          }
        }

        fs.open(cache(pid, "imiss_ihet").c_str(), ios::out | ios::binary);
        mpc.WriteToFile(imiss, fs);
        mpc.WriteToFile(ihet, fs);
        fs.close();

        cout << "Wrote results to cache" << endl;

        mpc.ProfilerPopState(false); // data_scan
      }

      mpc.ProfilerPushState("miss_filt");

      // Individual missingness filter
      cout << "Individual missing rate filter ... "; tic();
      ZZ_p imiss_ub = ZZ_p((long) (m1 * Param::IMISS_UB));
      mpc.LessThanPublic(ikeep, imiss, imiss_ub);
      cout << "done. "; toc();

      mpc.ProfilerPopState(true); // miss_filt
      mpc.ProfilerPushState("het_filt");

      // Individual heterozygosity filter
      cout << "Individual heterozygosity rate filter ... "; tic();
      ZZ_p ihet_ub_frac = DoubleToFP(Param::HET_UB, Param::NBIT_K, Param::NBIT_F);
      ZZ_p ihet_lb_frac = DoubleToFP(Param::HET_LB, Param::NBIT_K, Param::NBIT_F);

      // Number of observed SNPs per individual
      Vec<ZZ_p> m1_obs;
      Init(m1_obs, n0);
      if (pid > 0) {
        for (int i = 0; i < n0; i++) {
          m1_obs[i] = -imiss[i];
          if (pid == 1) {
            m1_obs[i] += m1;
          }
        }
      }

      Vec<ZZ_p> ihet_ub, ihet_lb;
      Init(ihet_ub, n0);
      Init(ihet_lb, n0);

      if (pid > 0) {
        for (int i = 0; i < n0; i++) {
          ihet_ub[i] = m1_obs[i] * ihet_ub_frac;
          ihet_lb[i] = m1_obs[i] * ihet_lb_frac;
          ihet[i] *= fp_one;
        }
      }

      Vec<ZZ_p> het_filt;
      mpc.LessThan(het_filt, ihet, ihet_ub);
      mpc.NotLessThan(tmp_vec, ihet, ihet_lb);
      mpc.MultElem(het_filt, het_filt, tmp_vec);

      mpc.MultElem(ikeep, ikeep, het_filt);
      het_filt.kill();
      cout << "done. "; toc();

      mpc.ProfilerPopState(true); // het_filt

      // Reveal samples to be filtered
      mpc.RevealSym(ikeep);

      if (pid == 2) {
        mpc.SendVec(ikeep, 0);
      } else if (pid == 0) {
        mpc.ReceiveVec(ikeep, 2, n0);
      }

      if (pid == 2) {
        ofs.open(outname("ikeep"));
        for (int i = 0; i < ikeep.length(); i++) {
          ofs << ikeep[i] << endl;
        }
        ofs.close();
      }
    }
  }

  mpc.ProfilerPopState(true); // ind_miss/het

  uint n1 = conv<uint>(Sum(ikeep));

  cout << "n1: " << n1 << ", " << "m1: " << m1 << endl;

  // Now calculate the "sub" n1 values for each dataset - ie number of individuals kept after filtering for each dataset
  vector<long> n1_vec;
  long rolling_n0 = 0;
  for (int i = 0; i < num_datasets; i++) {
    n1_vec.push_back(0);
    long sub_n0 = Param::NUM_INDS[i];
    for (int j = 0; j < sub_n0; j++) {
      if (ikeep[rolling_n0 + j] == 1) {
        n1_vec[i] = n1_vec[i] + 1;
      }
    }
    rolling_n0 += sub_n0;
  }

  cout << "Filtering phenotypes and covariates ... " << endl; tic();
  mpc.Filter(pheno, ikeep, n1);
  mpc.FilterRows(cov, ikeep, n1);
  cout << "done. "; toc();

  Vec<ZZ_p> ctrl;
  mpc.FlipBit(ctrl, pheno);

  Vec<ZZ_p> ctrl_mask;
  mpc.BeaverPartition(ctrl_mask, ctrl);

  Vec<ZZ_p> dosage_sum;
  Vec<ZZ_p> gmiss, gmiss_ctrl, dosage_sum_ctrl;
  Mat<ZZ_p> g_count_ctrl;
  ZZ_p n1_ctrl(0);

  Init(gmiss, m1);
  Init(gmiss_ctrl, m1);
  Init(dosage_sum, m1);
  Init(dosage_sum_ctrl, m1);
  Init(g_count_ctrl, 3, m1);

  mpc.ProfilerPushState("data_scan");

  if (exists(cache(pid, "geno_stats"))) {
    cout << "Genotype statistics cache found" << endl;

    ifs.open(cache(pid, "geno_stats").c_str(), ios::binary);
    mpc.ReadFromFile(gmiss, ifs, m1);
    mpc.ReadFromFile(gmiss_ctrl, ifs, m1);
    mpc.ReadFromFile(dosage_sum, ifs, m1);
    mpc.ReadFromFile(dosage_sum_ctrl, ifs, m1);
    mpc.ReadFromFile(g_count_ctrl, ifs, 3, m1);
    mpc.ReadFromFile(n1_ctrl, ifs);
    ifs.close();
  } else {
    cout << "Taking a pass to calculate genotype statistics:" << endl;

    // reduce batch size to avoid memory issues because of the overhead
    // in replicating dosage, g, miss, etc across all threads
    long bsize = Param::PITER_BATCH_SIZE / num_threads;
    long report_bsize = n1 / 10;

    // Loop over all datasets to calculate genotype statistics over all individuals
    #pragma omp parallel for num_threads(num_threads)
    for (int dataset_idx = 0; dataset_idx < num_datasets; dataset_idx++) {
      // Containers for batching the computation
      Vec< Mat<ZZ_p> > g, g_mask;
      Mat<ZZ_p> dosage, dosage_mask;
      Mat<ZZ_p> miss, miss_mask;
      Vec<ZZ_p> ctrl_vec, ctrl_mask_vec;
      g.SetLength(3);
      g_mask.SetLength(3);
      dosage.SetDims(bsize, m1);
      dosage_mask.SetDims(bsize, m1);
      miss.SetDims(bsize, m1);
      miss_mask.SetDims(bsize, m1);
      for (int k = 0; k < 3; k++) {
        g[k].SetDims(bsize, m1);
        g_mask[k].SetDims(bsize, m1);
      }
      ctrl_vec.SetLength(bsize);
      ctrl_mask_vec.SetLength(bsize);

      // Containers to store intermediate results for batching (costly) global updates
      Vec<ZZ_p> inner_gmiss, inner_gmiss_ctrl, inner_dosage_sum, inner_dosage_sum_ctrl;
      Mat<ZZ_p> inner_g_count_ctrl;
      ZZ_p inner_n1_ctrl(0);

      Init(inner_gmiss, m1);
      Init(inner_gmiss_ctrl, m1);
      Init(inner_dosage_sum, m1);
      Init(inner_dosage_sum_ctrl, m1);
      Init(inner_g_count_ctrl, 3, m1);

      tic();

      ifstream inner_ifs;
      inner_ifs.open(cache(pid, dataset_idx, "input_geno").c_str(), ios::binary);
      if (pid > 0) {
        mpc.ImportSeed(10, inner_ifs);
      } else {
        for (int p = 1; p <= 2; p++) {
          mpc.ImportSeed(10 + p, inner_ifs);
        }
      }

      // avoid error by re-setting the modulus within each thread
      ZZ base_p = conv<ZZ>(Param::BASE_P.c_str());
      ZZ_p::init(base_p);

      // Iterate over this dataset, considering only those individuals who have not been filtered out
      long offset = 0;
      long inner_ind = -1;
      for (int j = 0; j < dataset_idx; j++) {
        offset += n1_vec[j];
        inner_ind += Param::NUM_INDS[j];
      }
      long inner_n1 = n1_vec[dataset_idx];
      for (int i = 0; i < inner_n1; i++) {
        inner_ind++;

        mpc.ProfilerPushState("file_io/rng");

        Mat<ZZ_p> g0, g0_mask;
        Vec<ZZ_p> miss0, miss0_mask;

        while (ikeep[inner_ind] != 1) {
          if (pid > 0) {
            mpc.SkipData(inner_ifs, 3, m0); // g
            mpc.SkipData(inner_ifs, m0); // miss

            mpc.SwitchSeed(10);
            mpc.RandMat(g0_mask, 3, m0);
            mpc.RandVec(miss0_mask, m0);
            mpc.RestoreSeed();
          } else {
            for (int p = 1; p <= 2; p++) {
              mpc.SwitchSeed(10 + p);
              mpc.RandMat(g0_mask, 3, m0);
              mpc.RandVec(miss0_mask, m0);
              mpc.RestoreSeed();
            }
          }
          inner_ind++;
        }

        if (pid > 0) {
          mpc.ReadFromFile(g0, inner_ifs, 3, m0); // g
          mpc.ReadFromFile(miss0, inner_ifs, m0); // miss

          mpc.SwitchSeed(10);
          mpc.RandMat(g0_mask, 3, m0);
          mpc.RandVec(miss0_mask, m0);
          mpc.RestoreSeed();
        } else {
          Init(g0, 3, m0);
          Init(g0_mask, 3, m0);
          Init(miss0, m0);
          Init(miss0_mask, m0);
          Vec<ZZ_p> rand_vec;
          Mat<ZZ_p> rand_mat;

          for (int p = 1; p <= 2; p++) {
            mpc.SwitchSeed(10 + p);
            mpc.RandMat(rand_mat, 3, m0);
            mpc.RandVec(rand_vec, m0);
            mpc.RestoreSeed();

            g0_mask += rand_mat;
            miss0_mask += rand_vec;
          }
        }
        
        mpc.ProfilerPopState(false); // file_io/rng

        // Filter out loci that failed missing rate filter
        int inner_ind2 = 0;
        for (int j = 0; j < m0; j++) {
          if (gkeep1[j] == 1) {
            for (int k = 0; k < 3; k++) {
              g[k][i % bsize][inner_ind2] = g0[k][j];
              g_mask[k][i % bsize][inner_ind2] = g0_mask[k][j];
            }
            miss[i % bsize][inner_ind2] = miss0[j];
            miss_mask[i % bsize][inner_ind2] = miss0_mask[j];
            inner_ind2++;
          }
        }

        dosage[i % bsize] = g[1][i % bsize] + 2 * g[2][i % bsize];
        dosage_mask[i % bsize] = g_mask[1][i % bsize] + 2 * g_mask[2][i % bsize];

        ctrl_vec[i % bsize] = ctrl[i + offset];
        ctrl_mask_vec[i % bsize] = ctrl_mask[i + offset];

        // Update running sums
        if (pid > 0) {
          inner_n1_ctrl += ctrl_mask_vec[i % bsize];
          inner_gmiss += miss_mask[i % bsize];
          inner_dosage_sum += dosage_mask[i % bsize];

          if (pid == 1) {
              inner_n1_ctrl += ctrl_vec[i % bsize];
              inner_gmiss += miss[i % bsize];
              inner_dosage_sum += dosage[i % bsize];
          }
        }

        if (i % bsize == bsize - 1 || i == inner_n1 - 1) {
          if (i % bsize < bsize - 1) {
            int new_bsize = (i % bsize) + 1;
            for (int k = 0; k < 3; k++) {
              g[k].SetDims(new_bsize, m1);
              g_mask[k].SetDims(new_bsize, m1);
            }
            dosage.SetDims(new_bsize, m1);
            dosage_mask.SetDims(new_bsize, m1);
            miss.SetDims(new_bsize, m1);
            miss_mask.SetDims(new_bsize, m1);
            ctrl_vec.SetLength(new_bsize);
            ctrl_mask_vec.SetLength(new_bsize);
          }

          // Update running sums
          Vec<ZZ_p> tmp_gmiss_ctrl, tmp_dosage_sum_ctrl;
          Mat<ZZ_p> tmp_g_count_ctrl;
          Init(tmp_gmiss_ctrl, m1);
          Init(tmp_dosage_sum_ctrl, m1);
          Init(tmp_g_count_ctrl, 3, m1);

          mpc.BeaverMult(tmp_gmiss_ctrl, ctrl_vec, ctrl_mask_vec, miss, miss_mask);
          inner_gmiss_ctrl += tmp_gmiss_ctrl;

          mpc.BeaverMult(tmp_dosage_sum_ctrl, ctrl_vec, ctrl_mask_vec, dosage, dosage_mask);
          inner_dosage_sum_ctrl += tmp_dosage_sum_ctrl;

          for (int k = 0; k < 3; k++) {
            mpc.BeaverMult(tmp_g_count_ctrl[k], ctrl_vec, ctrl_mask_vec, g[k], g_mask[k]);
            inner_g_count_ctrl[k] += tmp_g_count_ctrl[k];
          }
        }

        if (i % report_bsize == 0 || i == inner_n1 - 1) {
          cout << "\t" << i << " / " << inner_n1 << ", "; toc(); tic();
        }
      }

      // Update global values - these operations must be atomic
      #pragma omp critical ( n1_ctrl_update )
        n1_ctrl += inner_n1_ctrl;
      #pragma omp critical ( gmiss_update )
        gmiss += inner_gmiss;
      #pragma omp critical ( dosage_sum_update )
        dosage_sum += inner_dosage_sum;
      #pragma omp critical ( gmiss_ctrl_update )
        gmiss_ctrl += inner_gmiss_ctrl;
      #pragma omp critical ( dosage_sum_ctrl_update )
        dosage_sum_ctrl += inner_dosage_sum_ctrl;
      #pragma omp critical ( g_count_ctrl_update )
        g_count_ctrl += inner_g_count_ctrl;

      inner_ifs.close();
    }

    mpc.BeaverReconstruct(gmiss_ctrl);
    mpc.BeaverReconstruct(dosage_sum_ctrl);
    mpc.BeaverReconstruct(g_count_ctrl);

    // Write to cache
    fs.open(cache(pid, "geno_stats").c_str(), ios::out | ios::binary);
    mpc.WriteToFile(gmiss, fs);
    mpc.WriteToFile(gmiss_ctrl, fs);
    mpc.WriteToFile(dosage_sum, fs);
    mpc.WriteToFile(dosage_sum_ctrl, fs);
    mpc.WriteToFile(g_count_ctrl, fs);
    mpc.WriteToFile(n1_ctrl, fs);
    fs.close();

    cout << "Wrote results to cache" << endl;
  }

  mpc.ProfilerPopState(true); // data_scan

  mpc.ProfilerPushState("maf/hwe");

  if (Param::DEBUG) {
    cout << "gmiss" << endl;
    mpc.Print(gmiss, 5);
    cout << "gmiss_ctrl" << endl;
    mpc.Print(gmiss_ctrl, 5);
    cout << "dosage_sum" << endl;
    mpc.Print(dosage_sum, 5);
    cout << "dosage_sum_ctrl" << endl;
    mpc.Print(dosage_sum_ctrl, 5);
    cout << "g_count_ctrl" << endl;
    for (int i = 0; i < 3; i++) {
      mpc.Print(g_count_ctrl[i], 5);
    }
  }

  mpc.ProfilerPushState("maf");

  // SNP MAF filter
  cout << "Locus minor allele frequency (MAF) filter ... " << endl; tic();
  ZZ_p maf_lb = DoubleToFP(Param::MAF_LB, Param::NBIT_K, Param::NBIT_F);
  ZZ_p maf_ub = DoubleToFP(Param::MAF_UB, Param::NBIT_K, Param::NBIT_F);

  Vec<ZZ_p> dosage_tot, dosage_tot_ctrl;
  if (pid > 0) {
    dosage_tot = -gmiss;
    dosage_tot_ctrl = -gmiss_ctrl;
    mpc.AddPublic(dosage_tot, ZZ_p(n1));
    mpc.Add(dosage_tot_ctrl, n1_ctrl);
    dosage_tot *= 2;
    dosage_tot_ctrl *= 2;
  } else {
    dosage_tot.SetLength(m1);
    dosage_tot_ctrl.SetLength(m1);
  }
  cout << "done. "; toc();

  cout << "Calculating MAFs ... " << endl; tic();
  Vec<ZZ_p> maf, maf_ctrl;
  if (exists(cache(pid, "maf"))) {
    cout << "maf cache found" << endl;
    ifs.open(cache(pid, "maf").c_str(), ios::binary);
    mpc.ReadFromFile(maf, ifs, dosage_tot.length());
    mpc.ReadFromFile(maf_ctrl, ifs, dosage_tot_ctrl.length());
    ifs.close();
  } else {
    mpc.ProfilerPushState("div");
    mpc.FPDivParallel(maf, dosage_sum, dosage_tot); 
    mpc.FPDivParallel(maf_ctrl, dosage_sum_ctrl, dosage_tot_ctrl); 
    mpc.ProfilerPopState(false); // div

    fs.open(cache(pid, "maf").c_str(), ios::out | ios::binary);
    mpc.WriteToFile(maf, fs);
    mpc.WriteToFile(maf_ctrl, fs);
    fs.close();
  }
  cout << "done. "; toc();

  Vec<ZZ_p> Maf, Maf_ctrl; // MAJOR allele freq
  if (pid > 0) {
    Maf = -maf;
    Maf_ctrl = -maf_ctrl;
    mpc.AddPublic(Maf, fp_one);
    mpc.AddPublic(Maf_ctrl, fp_one);
  } else {
    Maf.SetLength(m1);
    Maf_ctrl.SetLength(m1);
  }

  // Variance based on Bernoulli distribution over each allele
  Vec<ZZ_p> g_var_bern;
  mpc.MultElem(g_var_bern, maf, Maf);
  mpc.FastTrunc(g_var_bern);

  mpc.ProfilerPopState(true); // maf

  if (Param::DEBUG) {
    cout << "maf" << endl;
    mpc.PrintFP(maf, 5);
    cout << "maf_ctrl" << endl;
    mpc.PrintFP(maf_ctrl, 5);
  }

  Vec<ZZ_p> gkeep2;
  Init(gkeep2, m1);

  if (Param::SKIP_QC) {
    for (int i = 0; i < m1; i++) {
      gkeep2[i] = 1;
    }
    cout << "SNP MAF/HWE filters skipped" << endl;
  } else {
    bool history;
    if (pid == 2) {
      history = exists(outname("gkeep2"));
      mpc.SendBool(history, 0);
      mpc.SendBool(history, 1);
    } else {
      // ask P2 if gkeep2 has been computed before
      history = mpc.ReceiveBool(2);
    }

    if (history) {
      cout << "Using MAF/HWE filters from a previous run" << endl;
   
      if (pid == 2) {
        ifs.open(outname("gkeep2").c_str());
        for (int i = 0; i < m1; i++) {
          ifs >> gkeep2[i];
        }
        ifs.close();

        mpc.SendVec(gkeep2, 0);
        mpc.SendVec(gkeep2, 1);
      } else {
        mpc.ReceiveVec(gkeep2, 2, m1);
      }
    } else {
      
      mpc.ProfilerPushState("maf_filt");

      mpc.LessThanPublic(gkeep2, maf, maf_ub);
      mpc.NotLessThanPublic(tmp_vec, maf, maf_lb);
      mpc.MultElem(gkeep2, gkeep2, tmp_vec);

      mpc.ProfilerPopState(true); // maf_filt

      mpc.ProfilerPushState("hwe_filt");

      cout << "Locus Hardy-Weinberg equilibrium (HWE) filter ... " << endl; tic();
      ZZ_p hwe_ub = DoubleToFP(Param::HWE_UB, Param::NBIT_K, Param::NBIT_F); // p < 1e-7
      
      // Calculate expected genotype distribution in control group
      Mat<ZZ_p> g_exp_ctrl;
      Init(g_exp_ctrl, 3, m1);

      mpc.MultElem(g_exp_ctrl[0], Maf_ctrl, Maf_ctrl); // AA
      mpc.MultElem(g_exp_ctrl[1], Maf_ctrl, maf_ctrl); // Aa
      if (pid > 0) {
        g_exp_ctrl[1] *= 2;
      }
      mpc.MultElem(g_exp_ctrl[2], maf_ctrl, maf_ctrl); // aa

      for (int i = 0; i < 3; i++) {
        mpc.MultElem(g_exp_ctrl[i], g_exp_ctrl[i], dosage_tot_ctrl);
      }
      g_exp_ctrl *= twoinv; // dosage_tot_ctrl is twice the # individuals we actually want

      mpc.FastTrunc(g_exp_ctrl);

      cout << "\tCalculated expected genotype counts, "; toc(); tic();

      Vec<ZZ_p> hwe_chisq; 
      Init(hwe_chisq, m1);

      if (exists(cache(pid, "hwe"))) {
        cout << "HWE cache found" << endl;
        ifs.open(cache(pid, "hwe").c_str(), ios::binary);
        mpc.ReadFromFile(hwe_chisq, ifs, m1);
        ifs.close();
      } else {
        for (int i = 0; i < 3; i++) {
          Vec<ZZ_p> diff;
          if (pid > 0) {
            diff = fp_one * g_count_ctrl[i] - g_exp_ctrl[i];
          } else {
            diff.SetLength(m1);
          }

          mpc.MultElem(diff, diff, diff); // square
          mpc.FastTrunc(diff);

          mpc.ProfilerPushState("div");
          mpc.FPDivParallel(tmp_vec, diff, g_exp_ctrl[i]);
          mpc.ProfilerPopState(false); // div
          hwe_chisq += tmp_vec;

          cout << "\tChi-square test (" << i+1 << "/3), "; toc(); tic();
        }

        fs.open(cache(pid, "hwe").c_str(), ios::out | ios::binary);
        mpc.WriteToFile(hwe_chisq, fs);
        fs.close();
      }

      if (Param::DEBUG) {
        cout << "hwe" << endl;
        mpc.PrintFP(hwe_chisq, 5);
      }

      cout << "Applying hwe filter ... " << endl; tic();
      Vec<ZZ_p> hwe_filt;
      mpc.LessThanPublic(hwe_filt, hwe_chisq, hwe_ub);
      mpc.MultElem(gkeep2, gkeep2, hwe_filt);
      hwe_filt.kill();
      cout << "done. "; toc();

      // Reveal which SNPs to discard 
      mpc.RevealSym(gkeep2);
        
      if (pid == 2) {
        mpc.SendVec(gkeep2, 0);
      } else if (pid == 0) {
        mpc.ReceiveVec(gkeep2, 2, m1);
      }

      if (pid == 2) {
        ofs.open(outname("gkeep2").c_str());
        for (int i = 0; i < gkeep2.length(); i++) {
          ofs << gkeep2[i] << endl;
        }
        ofs.close();
      }

      mpc.ProfilerPopState(true); // hwe_filt
    }
  }

  uint m2 = conv<uint>(Sum(gkeep2));
  cout << "n1: " << n1 << ", " << "m2: " << m2 << endl;

  cout << "Filtering genotype statistics" << endl; tic();
  mpc.Filter(g_var_bern, gkeep2, m2);
  mpc.Filter(maf, gkeep2, m2);
  FilterVec(snp_pos, gkeep2);
  cout << "done. "; toc();

  gmiss.kill();
  gmiss_ctrl.kill();
  dosage_sum.kill();
  dosage_sum_ctrl.kill();
  g_count_ctrl.kill();

  mpc.ProfilerPopState(false); // maf/hwe
  mpc.ProfilerPopState(true); // qc
  mpc.ProfilerPushState("std_param");

  Vec<ZZ_p> g_std_bern_inv;
  if (exists(cache(pid, "stdinv_bern"))) {
    cout << "Genotype standard deviation cache found" << endl;

    ifs.open(cache(pid, "stdinv_bern").c_str(), ios::binary);
    mpc.ReadFromFile(g_std_bern_inv, ifs, g_var_bern.length());
    ifs.close();

  } else {
    cout << "Calculating genotype standard deviations (inverse)" << endl; tic();

    mpc.ProfilerPushState("sqrt");
    mpc.FPSqrtParallel(tmp_vec, g_std_bern_inv, g_var_bern);
    mpc.ProfilerPopState(false); // sqrt

    cout << "done. "; toc();

    fs.open(cache(pid, "stdinv_bern").c_str(), ios::out | ios::binary);
    mpc.WriteToFile(g_std_bern_inv, fs);
    fs.close();
  }

  if (Param::DEBUG) {
    cout << "g_std_bern_inv" << endl;
    mpc.PrintFP(g_std_bern_inv, 5);
  }

  Vec<ZZ_p> g_mean;
  if (pid > 0) {
    g_mean = 2 * maf;
  } else {
    g_mean.SetLength(m2);
  }

  mpc.ProfilerPopState(true); // std_param

  cout << "Starting population stratification analysis" << endl;

  mpc.ProfilerPushState("pop_strat");
  mpc.ProfilerPushState("select_snp");

  Vec<int8_t> selected; // 1 selected, 0 unselected, -1 TBD
  Vec<bool> to_process;
  selected.SetLength(m2);
  to_process.SetLength(m2);

  for (int i = 0; i < m2; i++) {
    selected[i] = -1;
  }

  ZZ dist_thres(Param::LD_DIST_THRES);
 
  ZZ prev(-1);
  for (int i = 0; i < m2; i++) {
    selected[i] = 0;
    if (prev < 0 || snp_pos[i] - prev > dist_thres) {
      selected[i] = 1;
      prev = snp_pos[i];
    }
  }

  // At this point "selected" contains the SNP filter for PCA, shared across all parties
  uint32_t m3 = 0;
  for (int i = 0; i < selected.length(); i++) {
    if (selected[i] == 1) {
      m3++;
    }
  }

  cout << "SNP selection complete: " << m3 << " / " << m2 << " selected" << endl;
  mpc.ProfilerPopState(false); // select_snp
  mpc.ProfilerPushState("reduce_file");

  // Cache the reduced G for PCA
  if (exists(cache(pid, "pca_input"))) {
    cout << "pca_input cache found" << endl;
  } else {
    Vec<bool> gkeep3;
    gkeep3.SetLength(m0);
    for (int j = 0; j < m0; j++) {
      gkeep3[j] = gkeep1[j] == 1;
    }

    ind = 0;
    for (int j = 0; j < m0; j++) {
      if (gkeep3[j]) {
        gkeep3[j] = gkeep2[ind] == 1;
        ind++;
      }
    }

    ind = 0;
    for (int j = 0; j < m0; j++) {
      if (gkeep3[j]) {
        gkeep3[j] = selected[ind] == 1;
        ind++;
      }
    }

    long bsize = n1 / 10;

    cout << "Caching input data for PCA:" << endl;

    // Loop over all datasets
    #pragma omp parallel for num_threads(num_threads)
    for (int dataset_idx = 0; dataset_idx < num_datasets; dataset_idx++) {
      tic();
      ifstream inner_ifs;
      inner_ifs.open(cache(pid, dataset_idx, "input_geno").c_str(), ios::binary);
      if (pid > 0) {
        mpc.ImportSeed(10, inner_ifs);
      } else {
        for (int p = 1; p <= 2; p++) {
          mpc.ImportSeed(10 + p, inner_ifs);
        }
      }

      // Improve performance by writing the pca data for each chunk to a separate cache file in parallel
      fstream inner_fs;
      inner_fs.open(cache(pid, dataset_idx, "pca_input").c_str(), ios::out | ios::binary);

      // avoid error by re-setting the modulus within each thread
      ZZ base_p = conv<ZZ>(Param::BASE_P.c_str());
      ZZ_p::init(base_p);

      // Iterate over this dataset, considering only those individuals who have not been filtered out
      long offset = 0;
      long inner_ind = -1;
      for (int j = 0; j < dataset_idx; j++) {
        offset += n1_vec[j];
        inner_ind += Param::NUM_INDS[j];
      }
      long inner_n1 = n1_vec[dataset_idx];
      for (int i = offset; i < offset + inner_n1; i++) {
        inner_ind++;

        mpc.ProfilerPushState("file_io/rng");

        Mat<ZZ_p> g0, g0_mask;
        Vec<ZZ_p> miss0, miss0_mask;

        while (ikeep[inner_ind] != 1) {
          if (pid > 0) {
            mpc.SkipData(inner_ifs, 3, m0); // g
            mpc.SkipData(inner_ifs, m0); // miss

            mpc.SwitchSeed(10);
            mpc.RandMat(g0_mask, 3, m0);
            mpc.RandVec(miss0_mask, m0);
            mpc.RestoreSeed();
          } else {
            for (int p = 1; p <= 2; p++) {
              mpc.SwitchSeed(10 + p);
              mpc.RandMat(g0_mask, 3, m0);
              mpc.RandVec(miss0_mask, m0);
              mpc.RestoreSeed();
            }
          }
          inner_ind++;
        }

        if (pid > 0) {
          mpc.ReadFromFile(g0, inner_ifs, 3, m0); // g
          mpc.ReadFromFile(miss0, inner_ifs, m0); // miss

          mpc.SwitchSeed(10);
          mpc.RandMat(g0_mask, 3, m0);
          mpc.RandVec(miss0_mask, m0);
          mpc.RestoreSeed();
        } else {
          Init(g0, 3, m0);
          Init(g0_mask, 3, m0);
          Init(miss0, m0);
          Init(miss0_mask, m0);
          Vec<ZZ_p> rand_vec;
          Mat<ZZ_p> rand_mat;

          for (int p = 1; p <= 2; p++) {
            mpc.SwitchSeed(10 + p);
            mpc.RandMat(rand_mat, 3, m0);
            mpc.RandVec(rand_vec, m0);
            mpc.RestoreSeed();

            g0_mask += rand_mat;
            miss0_mask += rand_vec;
          }
        }
        
        mpc.ProfilerPopState(false); // file_io/rng

        // Filter out loci that failed missing rate filter
        Mat<ZZ_p> g, g_mask;
        Vec<ZZ_p> miss, miss_mask;
        g.SetDims(3, m3);
        g_mask.SetDims(3, m3);
        miss.SetLength(m3);
        miss_mask.SetLength(m3);

        int ind2 = 0;
        for (int j = 0; j < m0; j++) {
          if (gkeep3[j]) {
            for (int k = 0; k < 3; k++) {
              g[k][ind2] = g0[k][j];
              g_mask[k][ind2] = g0_mask[k][j];
            }
            miss[ind2] = miss0[j];
            miss_mask[ind2] = miss0_mask[j];
            ind2++;
          }
        }

        Vec<ZZ_p> dosage, dosage_mask;
        dosage = g[1] + 2 * g[2];
        dosage_mask = g_mask[1] + 2 * g_mask[2];

        mpc.BeaverWriteToFile(dosage, dosage_mask, inner_fs);
        mpc.BeaverWriteToFile(miss, miss_mask, inner_fs);

        if ((i - offset + 1) % bsize == 0 || (i - offset) == inner_n1 - 1) {
          cout << "\t" << i+1 << " / " << n1 << ", "; toc(); tic();
        }
      }
      inner_ifs.close();
      inner_fs.close();
    }
  }

  mpc.ProfilerPopState(false); // reduce_file

  Vec<ZZ_p> g_mean_pca = g_mean;
  mpc.Filter(g_mean_pca, selected, m3);

  Vec<ZZ_p> g_stdinv_pca = g_std_bern_inv;
  mpc.Filter(g_stdinv_pca, selected, m3);

  Vec<ZZ_p> g_mean_pca_mask, g_stdinv_pca_mask;
  mpc.BeaverPartition(g_mean_pca_mask, g_mean_pca);
  mpc.BeaverPartition(g_stdinv_pca_mask, g_stdinv_pca);

  /* Pass 2: Random sketch */
  Mat<ZZ_p> Y_cur;
  Init(Y_cur, kp, m3);

  if (exists(cache(pid, "sketch"))) {

    cout << "sketch cache found" << endl;
    ifs.open(cache(pid, "sketch").c_str(), ios::in | ios::binary);
    ifs >> kp;
    mpc.ReadFromFile(Y_cur, ifs, kp, m3);
    ifs.close();

  } else {

    mpc.ProfilerPushState("rand_proj");

    cout << "Random Projection ..." << endl; tic();
    
    Mat<ZZ_p> Y_cur_adj;
    Init(Y_cur_adj, kp, m3);

    Vec<int> bucket_count;
    bucket_count.SetLength(kp);
    for (int i = 0; i < kp; i++) {
      bucket_count[i] = 0;
    }

    #pragma omp parallel for num_threads(num_threads)
    for (int dataset_idx = 0; dataset_idx < num_datasets; dataset_idx++) {
      ifstream inner_ifs;
      inner_ifs.open(cache(pid, dataset_idx, "pca_input").c_str(), ios::in | ios::binary);

      // Containers to store intermediate results for batching (costly) global updates
      Mat<ZZ_p> inner_Y_cur, inner_Y_cur_adj;
      Init(inner_Y_cur, kp, m3);
      Init(inner_Y_cur_adj, kp, m3);

      Vec<int> inner_bucket_count;
      inner_bucket_count.SetLength(kp);
      for (int i = 0; i < kp; i++) {
        inner_bucket_count[i] = 0;
      }

      long inner_n1 = n1_vec[dataset_idx];
      for (int cur = 0; cur < inner_n1; cur++) {
        // Count sketch (use global PRG)
        mpc.SwitchSeed(-1);
        long bucket_index = mpc.RandBnd(kp);
        long rand_sign = mpc.RandBnd(2) * 2 - 1;
        mpc.RestoreSeed();

        Vec<ZZ_p> g, g_mask, miss, miss_mask;
        mpc.BeaverReadFromFile(g, g_mask, inner_ifs, m3);
        mpc.BeaverReadFromFile(miss, miss_mask, inner_ifs, m3);

        // Flip miss bits so it points to places where g_mean should be subtracted
        mpc.BeaverFlipBit(miss, miss_mask);

        // Update running sum
        if (pid > 0) {
          inner_Y_cur[bucket_index] += rand_sign * g_mask;
          if (pid == 1) {
            inner_Y_cur[bucket_index] += rand_sign * g;
          }
        }

        // Update adjustment factor
        miss *= rand_sign;
        miss_mask *= rand_sign;
        mpc.BeaverMultElem(inner_Y_cur_adj[bucket_index], miss, miss_mask, g_mean_pca, g_mean_pca_mask);

        inner_bucket_count[bucket_index]++;
      }
      inner_ifs.close();

       // Update global values - these operations must be atomic
      #pragma omp critical ( Y_cur_update )
        Y_cur += inner_Y_cur;
      #pragma omp critical ( Y_cur_adj_update )
        Y_cur_adj += inner_Y_cur_adj;
      #pragma omp critical ( bucket_count_update )
        for (int i = 0; i < kp; i++) {
          bucket_count[i] = bucket_count[i] + inner_bucket_count[i];
        }
    }

    // Subtract the adjustment factor
    mpc.BeaverReconstruct(Y_cur_adj);
    if (pid > 0) {
      Y_cur = fp_one * Y_cur - Y_cur_adj;
    }
    Y_cur_adj.kill();

    cout << "done. "; toc();

    if (Param::DEBUG) {
      cout << "Y_cur" << endl;
      mpc.PrintFP(Y_cur[0], 5);
      cout << "g_mean_pca" << endl;
      mpc.PrintBeaverFP(g_mean_pca, g_mean_pca_mask, 10);
      cout << "g_stdinv_pca" << endl;
      mpc.PrintBeaverFP(g_stdinv_pca, g_stdinv_pca_mask, 10);
    }

    // Get rid of empty buckets and normalize nonempty ones
    int empty_slot = 0;
    for (int i = 0; i < kp; i++) {
      if (bucket_count[i] > 0) {
        ZZ_p fp_count_inv = DoubleToFP(1 / ((double) bucket_count[i]), Param::NBIT_K, Param::NBIT_F);
        Y_cur[empty_slot] = Y_cur[i] * fp_count_inv;
        empty_slot++;
      }
    }
    kp = empty_slot;
    cout << "kp: " << kp << endl;
    Y_cur.SetDims(kp, m3);
    mpc.FastTrunc(Y_cur);

    mpc.ProfilerPopState(true); // rand_proj

    fs.open(cache(pid, "sketch").c_str(), ios::out | ios::binary);
    fs << kp;
    if (pid > 0) {
      mpc.WriteToFile(Y_cur, fs);
    }
    fs.close();
  }

  mpc.ProfilerPushState("power_iter");

  Mat<ZZ_p> Y_cur_mask;
  mpc.BeaverPartition(Y_cur_mask, Y_cur);

  cout << "Initial sketch obtained, starting power iteration (num iter = " << Param::NUM_POWER_ITER << ")" << endl;
  tic();

  Mat<ZZ_p> gQ;

  if (exists(cache(pid, "piter"))) {

    cout << "piter cache found" << endl;
    ifs.open(cache(pid, "piter").c_str(), ios::in | ios::binary);
    mpc.ReadFromFile(gQ, ifs, n1, kp);
    ifs.close();
    
  } else {

    // Divide by standard deviation
    Mat<ZZ_p> Y;
    Init(Y, kp, m3);

    for (int i = 0; i < kp; i++) {
      mpc.BeaverMultElem(Y[i], Y_cur[i], Y_cur_mask[i], g_stdinv_pca, g_stdinv_pca_mask);
    }
    Y_cur.kill();
    Y_cur_mask.kill();

    mpc.BeaverReconstruct(Y);
    mpc.FastTrunc(Y);

    /* Calculate orthonormal bases of Y */
    cout << "Initial orthonormal basis ... "; tic();
    Mat<ZZ_p> Q;
    mpc.ProfilerPushState("qr_m");
    mpc.OrthonormalBasis(Q, Y);
    mpc.ProfilerPopState(false); // qr_m
    Y.kill();
    cout << "done. "; toc();

    Mat<ZZ_p> gQ_adj;
    Mat<ZZ_p> Q_mask;
    Mat<ZZ_p> Q_scaled, Q_scaled_mask;
    Mat<ZZ_p> Q_scaled_gmean, Q_scaled_gmean_mask;

    /* Power iteration */
    for (int pit = 0; pit <= Param::NUM_POWER_ITER; pit++) {
      /* This section is ran before each iteration AND once after all iterations */
      mpc.BeaverPartition(Q_mask, Q);

      // Normalize Q by standard deviations
      Init(Q_scaled, kp, m3);
      for (int i = 0; i < kp; i++) {
        mpc.BeaverMultElem(Q_scaled[i], Q[i], Q_mask[i], g_stdinv_pca, g_stdinv_pca_mask);
      }
      mpc.BeaverReconstruct(Q_scaled);
      mpc.FastTrunc(Q_scaled);

      mpc.BeaverPartition(Q_scaled_mask, Q_scaled);

      // Pre-multiply with g_mean to simplify calculation of centering matrix
      Init(Q_scaled_gmean, kp, m3);
      for (int i = 0; i < kp; i++) {
        mpc.BeaverMultElem(Q_scaled_gmean[i], Q_scaled[i], Q_scaled_mask[i],
                           g_mean_pca, g_mean_pca_mask);
      }
      mpc.BeaverReconstruct(Q_scaled_gmean);
      mpc.FastTrunc(Q_scaled_gmean);

      mpc.Transpose(Q_scaled); // m3-by-kp
      transpose(Q_scaled_mask, Q_scaled_mask); // m3-by-kp, unlike mpc.Transpose, P0 also transposes
      mpc.Transpose(Q_scaled_gmean); // m3-by-kp
      mpc.BeaverPartition(Q_scaled_gmean_mask, Q_scaled_gmean);

      // reduce batch size to avoid memory issues when replicating variables across threads
      long bsize = Param::PITER_BATCH_SIZE / num_threads;
      
      /* Pass 1 */
      Init(gQ, n1, kp);
      Init(gQ_adj, n1, kp);

      mpc.ProfilerPushState("data_scan0");

      if (pit == 0) {
        cout << "Iter 1: Data Scan 1 ... "; tic();
      }

      #pragma omp parallel for num_threads(num_threads)
      for (int dataset_idx = 0; dataset_idx < num_datasets; dataset_idx++) {
        mpc.ProfilerPushState("file_io");
        ifstream inner_ifs;
        inner_ifs.open(cache(pid, dataset_idx, "pca_input").c_str(), ios::in | ios::binary);

        Mat<ZZ_p> g, g_mask, miss, miss_mask;
        Init(g, bsize, m3);
        Init(g_mask, bsize, m3);
        Init(miss, bsize, m3);
        Init(miss_mask, bsize, m3);

        long offset = 0;
        for (int j = 0; j < dataset_idx; j++) {
          offset += n1_vec[j];
        }
        long inner_n1 = n1_vec[dataset_idx];
        for (int cur = 0; cur < inner_n1; cur++) {
          mpc.BeaverReadFromFile(g[cur % bsize], g_mask[cur % bsize], inner_ifs, m3);
          mpc.BeaverReadFromFile(miss[cur % bsize], miss_mask[cur % bsize], inner_ifs, m3);
          mpc.BeaverFlipBit(miss[cur % bsize], miss_mask[cur % bsize]);

          if (cur % bsize == bsize - 1 || cur == inner_n1 - 1) {
            mpc.ProfilerPopState(false); // file_io
            int new_bsize = bsize;
            if (cur % bsize < bsize - 1) {
              new_bsize = (cur % bsize) + 1;
              g.SetDims(new_bsize, m3);
              g_mask.SetDims(new_bsize, m3);
              miss.SetDims(new_bsize, m3);
              miss_mask.SetDims(new_bsize, m3);
            }

            Mat<ZZ_p> result;
            Init(result, new_bsize, kp);
            mpc.BeaverMult(result, g, g_mask, Q_scaled, Q_scaled_mask);
            for (int i = 0; i < new_bsize; i++) {
              gQ[(cur+offset)-(new_bsize-1)+i] = result[i];
            }

            Init(result, new_bsize, kp);
            mpc.BeaverMult(result, miss, miss_mask, Q_scaled_gmean, Q_scaled_gmean_mask);
            for (int i = 0; i < new_bsize; i++) {
              gQ_adj[(cur+offset)-(new_bsize-1)+i] = result[i];
            }

            if (cur < inner_n1 - 1) mpc.ProfilerPushState("file_io");
          }
        }
        ifs.close();
        
        g.kill();
        g_mask.kill();
        miss.kill();
        miss_mask.kill();
      }
      if (pit == 0) {
        cout << "done. "; toc();
      }

      mpc.BeaverReconstruct(gQ);
      mpc.BeaverReconstruct(gQ_adj);
      if (pid > 0) {
        gQ -= gQ_adj;
      }

      mpc.ProfilerPopState(false); // data_scan1

      if (pit == Param::NUM_POWER_ITER) { // Quit if all iterations are performed
        break;
      }

      if (pit == 0) {
        cout << "Iter 1: OrthonormalBasis 1 ... "; tic();
      }
      mpc.Transpose(gQ); // kp-by-n1
      mpc.ProfilerPushState("qr_n");
      mpc.OrthonormalBasis(Q, gQ);
      mpc.ProfilerPopState(false); // qr_n
      mpc.Transpose(Q); // n1-by-kp
      if (pit == 0) {
        cout << "done. "; toc();
      }

      mpc.BeaverPartition(Q_mask, Q);

      Init(gQ, kp, m3);
      Init(gQ_adj, kp, m3);

      mpc.ProfilerPushState("data_scan2");

      // Pass 2
      if (pit == 0) {
        cout << "Iter 1: Data Scan 2 ... "; tic();
      }

      #pragma omp parallel for num_threads(num_threads)
      for (int dataset_idx = 0; dataset_idx < num_datasets; dataset_idx++) {
        mpc.ProfilerPushState("file_io");
        ifstream inner_ifs;
        inner_ifs.open(cache(pid, dataset_idx, "pca_input").c_str(), ios::in | ios::binary);

        Mat<ZZ_p> g, g_mask, miss, miss_mask, Qsub, Qsub_mask;;
        Init(g, bsize, m3);
        Init(g_mask, bsize, m3);
        Init(miss, bsize, m3);
        Init(miss_mask, bsize, m3);
        Init(Qsub, bsize, kp);
        Init(Qsub_mask, bsize, kp);

        // Containers to store intermediate results for batching (costly) global updates
        Mat<ZZ_p> inner_gQ, inner_gQ_adj;
        Init(inner_gQ, kp, m3);
        Init(inner_gQ_adj, kp, m3);

        long offset = 0;
        for (int j = 0; j < dataset_idx; j++) {
          offset += n1_vec[j];
        }
        long inner_n1 = n1_vec[dataset_idx];
        for (int cur = 0; cur < inner_n1; cur++) {
          mpc.BeaverReadFromFile(g[cur % bsize], g_mask[cur % bsize], inner_ifs, m3);
          mpc.BeaverReadFromFile(miss[cur % bsize], miss_mask[cur % bsize], inner_ifs, m3);
          mpc.BeaverFlipBit(miss[cur % bsize], miss_mask[cur % bsize]);

          Qsub[cur % bsize] = Q[cur + offset];
          Qsub_mask[cur % bsize] = Q_mask[cur + offset];

          if (cur % bsize == bsize - 1 || cur == inner_n1 - 1) {
            mpc.ProfilerPopState(false); // file_io
            int new_bsize = bsize;
            if (cur % bsize < bsize - 1) {
              new_bsize = (cur % bsize) + 1;
              g.SetDims(new_bsize, m3);
              g_mask.SetDims(new_bsize, m3);
              miss.SetDims(new_bsize, m3);
              miss_mask.SetDims(new_bsize, m3);
              Qsub.SetDims(new_bsize, kp);
              Qsub_mask.SetDims(new_bsize, kp);
            }

            mpc.Transpose(Qsub);
            transpose(Qsub_mask, Qsub_mask);

            mpc.BeaverMult(inner_gQ, Qsub, Qsub_mask, g, g_mask);
            mpc.BeaverMult(inner_gQ_adj, Qsub, Qsub_mask, miss, miss_mask);

            Qsub.SetDims(bsize, kp);
            Qsub_mask.SetDims(bsize, kp);

            if (cur < inner_n1 - 1) mpc.ProfilerPushState("file_io");
          }
        }
        ifs.close();

        // Update global values - these operations must be atomic
        #pragma omp critical ( gQ_update )
          gQ += inner_gQ;
        #pragma omp critical ( gQ_adj_update )
          gQ_adj += inner_gQ_adj;
        
        g.kill();
        g_mask.kill();
        miss.kill();
        miss_mask.kill();
        Qsub.kill();
        Qsub_mask.kill();
        inner_gQ.kill();
        inner_gQ_adj.kill();
      }
      if (pit == 0) {
        cout << "done. "; toc();
      }

      mpc.BeaverReconstruct(gQ);
      mpc.BeaverReconstruct(gQ_adj);

      mpc.ProfilerPopState(false); // data_scan2

      Mat<ZZ_p> gQ_adj_mask;
      mpc.BeaverPartition(gQ_adj_mask, gQ_adj);

      Mat<ZZ_p> gQ_adj_gmean;
      Init(gQ_adj_gmean, kp, m3);
      for (int i = 0; i < kp; i++) {
        mpc.BeaverMultElem(gQ_adj_gmean[i], gQ_adj[i], gQ_adj_mask[i],
                           g_mean_pca, g_mean_pca_mask);
      }
      mpc.BeaverReconstruct(gQ_adj_gmean);
      mpc.FastTrunc(gQ_adj_gmean);

      if (pid > 0) {
        gQ -= gQ_adj_gmean;
      }
      gQ_adj_gmean.kill();

      Mat<ZZ_p> gQ_mask;
      mpc.BeaverPartition(gQ_mask, gQ);

      Mat<ZZ_p> gQ_scaled;
      gQ_scaled.SetDims(kp, m3);
      clear(gQ_scaled);
      for (int i = 0; i < kp; i++) {
        mpc.BeaverMultElem(gQ_scaled[i], gQ[i], gQ_mask[i], g_stdinv_pca, g_stdinv_pca_mask);
      }
      mpc.BeaverReconstruct(gQ_scaled);
      mpc.FastTrunc(gQ_scaled);

      mpc.ProfilerPushState("qr_m");
      if (pit == 0) {
        cout << "Iter 1: Orthonormal Basis 2 ... "; tic();
      }
      mpc.OrthonormalBasis(Q, gQ_scaled);
      if (pit == 0) {
        cout << "done. "; toc();
      }
      mpc.ProfilerPopState(false); // qr_m

      if (pit != 0) {
        cout << "Iter " << pit + 1 << " complete, "; toc();
      }
      tic();
    }

    fs.open(cache(pid, "piter").c_str(), ios::out | ios::binary);
    if (pid > 0) {
      mpc.WriteToFile(gQ, fs);
    }
    fs.close();
  }

  mpc.ProfilerPopState(true); // power_iter
  cout << "Power iteration complete" << endl;

  Mat<ZZ_p> Z = gQ;
  gQ.kill();

  cout << "Data projected to subspace" << endl;
  if (Param::DEBUG) {
    cout << "Z" << endl;
    mpc.PrintFP(Z[0], 5);
  }

  Mat<ZZ_p> V;
  Init(V, k, n1);

  /* Eigendecomposition */
  if (exists(cache(pid, "eigen"))) {

    cout << "eigen cache found" << endl;
    ifs.open(cache(pid, "eigen").c_str(), ios::binary);
    mpc.ReadFromFile(V, ifs, k, n1);
    ifs.close();

  } else {

    ZZ_p fp_m2_inv = DoubleToFP(1 / ((double) m2), Param::NBIT_K, Param::NBIT_F);
    Z *= fp_m2_inv;
    mpc.FastTrunc(Z);

    mpc.Transpose(Z); // kp-by-n1

    Mat<ZZ_p> Z_mask;
    mpc.BeaverPartition(Z_mask, Z);

    /* Form covariance matrix */
    Mat<ZZ_p> Z_gram;
    Init(Z_gram, kp, kp);
    for (int i = 0; i < kp; i++) {
      mpc.BeaverMult(Z_gram[i], Z, Z_mask, Z[i], Z_mask[i]);
    }
    mpc.BeaverReconstruct(Z_gram);
    mpc.FastTrunc(Z_gram);

    cout << "Constructed reduced eigenvalue problem" << endl;

    if (Param::DEBUG) {
      cout << "Z_gram" << endl;
      for (int i = 0; i < 3; i++) {
        mpc.PrintFP(Z_gram[i], 5);
      }
    }

    mpc.ProfilerPushState("eigen_solve");

    Mat<ZZ_p> U;
    Vec<ZZ_p> L;
    cout << "Eigenvector decomposition ... " << endl; tic();
    mpc.EigenDecomp(U, L, Z_gram);
    cout << "done. "; toc();
    Z_gram.kill();

    // Select top eigenvectors and eigenvalues
    U.SetDims(k, kp);
    L.SetLength(k);

    cout << "Selected K eigenvectors" << endl;
    mpc.ProfilerPopState(false); // eigen_solve

    if (Param::DEBUG) {
      for (int i = 0; i < 3; i++) {
        mpc.PrintFP(U[i], 5);
      }
      cout << "Eigenvalues" << endl;
      mpc.PrintFP(L, 5);
    }

    // Recover singular vectors
    Mat<ZZ_p> U_mask;
    mpc.BeaverPartition(U_mask, U);

    mpc.BeaverMultMat(V, U, U_mask, Z, Z_mask);
    U.kill();
    U_mask.kill();
    Z_mask.kill();
    mpc.BeaverReconstruct(V);
    mpc.FastTrunc(V);

    fs.open(cache(pid, "eigen").c_str(), ios::out | ios::binary);
    if (pid > 0) {
      mpc.WriteToFile(V, fs);
    }
    fs.close();

  }

  Z.kill();

  mpc.ProfilerPopState(true); // pop_strat

  mpc.ProfilerPushState("assoc_test");
  mpc.ProfilerPushState("covar");

  // Concatenate covariate matrix and jointly orthogonalize
  mpc.Transpose(cov);
  V.SetDims(k + Param::NUM_COVS, n1);
  if (pid > 0) {
    for (int i = 0; i < Param::NUM_COVS; i++) {
      V[k + i] = cov[i] * fp_one;
    }
  }
  cov.kill();
  mpc.OrthonormalBasis(V, V);

  Mat<ZZ_p> V_mask;
  mpc.BeaverPartition(V_mask, V);

  cout << "Bases for top singular vectors and covariates calculated" << endl;
  mpc.ProfilerPopState(false); // covar

  if (Param::DEBUG) {
    mpc.PrintBeaverFP(V[0], V_mask[0], 5);
  }

  /* Pass 4: Calculate GWAS statistics */

  Vec<ZZ_p> pheno_mask;
  mpc.BeaverPartition(pheno_mask, pheno);

  Vec<ZZ_p> Vp;
  Init(Vp, k + Param::NUM_COVS);
  mpc.BeaverMult(Vp, V, V_mask, pheno, pheno_mask);
  mpc.BeaverReconstruct(Vp);

  Vec<ZZ_p> Vp_mask;
  mpc.BeaverPartition(Vp_mask, Vp);
  
  Vec<ZZ_p> VVp;
  Init(VVp, n1);
  mpc.BeaverMult(VVp, Vp, Vp_mask, V, V_mask);
  mpc.BeaverReconstruct(VVp);
  mpc.FastTrunc(VVp);

  Vec<ZZ_p> VVp_mask;
  mpc.BeaverPartition(VVp_mask, VVp);

  Vec<ZZ_p> p_hat, p_hat_mask;
  p_hat = fp_one * pheno - VVp;
  p_hat_mask = fp_one * pheno_mask - VVp_mask;

  Vp.kill();
  Vp_mask.kill();
  VVp.kill();
  VVp_mask.kill();

  cout << "Phenotypes corrected" << endl;

  Vec<ZZ_p> V_sum, V_sum_mask;
  Init(V_sum, k + Param::NUM_COVS);
  Init(V_sum_mask, k + Param::NUM_COVS);
  for (int i = 0; i < k + Param::NUM_COVS; i++) {
    for (int j = 0; j < n1; j++) {
      V_sum[i] += V[i][j];
      V_sum_mask[i] += V_mask[i][j];
    }
  }

  Vec<ZZ_p> u;
  Init(u, n1);
  mpc.BeaverMult(u, V_sum, V_sum_mask, V, V_mask);
  mpc.BeaverReconstruct(u);
  mpc.FastTrunc(u);
  if (pid > 0) {
    u *= -1;
    mpc.AddPublic(u, fp_one);
  }

  Vec<ZZ_p> u_mask;
  mpc.BeaverPartition(u_mask, u);

  if (Param::DEBUG) {
    cout << "u" << endl;
    mpc.PrintBeaverFP(u, u_mask, 10);
  }

  cout << "Allocating sx, sxx, sxp, B ... ";

  Vec<ZZ_p> sx, sxx, sxp;
  Mat<ZZ_p> B;
  Init(sx, m2);
  Init(sxx, m2);
  Init(sxp, m2);
  Init(B, k + Param::NUM_COVS, m2);

  cout << "done.";

  mpc.ProfilerPushState("data_scan");

  if (exists(cache(pid, "gwas_stats"))) {
    cout << "GWAS statistics cache found" << endl;
    ifs.open(cache(pid, "gwas_stats").c_str(), ios::binary);
    mpc.ReadFromFile(sx, ifs, m2);
    mpc.ReadFromFile(sxx, ifs, m2);
    mpc.ReadFromFile(sxp, ifs, m2);
    mpc.ReadFromFile(B, ifs, k + Param::NUM_COVS, m2);
    ifs.close();

  } else {

    long bsize = Param::PITER_BATCH_SIZE / num_threads;

    mpc.Transpose(V); // n1-by-(k + NUM_COVS)
    transpose(V_mask, V_mask);

    Vec<bool> gkeep3;
    gkeep3.SetLength(m0);
    for (int j = 0; j < m0; j++) {
      gkeep3[j] = gkeep1[j] == 1;
    }

    ind = 0;
    for (int j = 0; j < m0; j++) {
      if (gkeep3[j]) {
        gkeep3[j] = gkeep2[ind] == 1;
        ind++;
      }
    }

    tic();
    mpc.ProfilerPushState("file_io/rng");
    cout << "GWAS pass:" << endl;
    
    #pragma omp parallel for num_threads(num_threads)
    for (int dataset_idx = 0; dataset_idx < Param::NUM_INDS.size(); dataset_idx++) {
      Mat<ZZ_p> dosage, dosage_mask;
      Init(dosage, bsize, m2);
      Init(dosage_mask, bsize, m2);

      Vec<ZZ_p> u_vec, u_mask_vec, p_hat_vec, p_hat_mask_vec;
      Init(u_vec, bsize);
      Init(u_mask_vec, bsize);
      Init(p_hat_vec, bsize);
      Init(p_hat_mask_vec, bsize);

      Mat<ZZ_p> V_sub, V_mask_sub;
      Init(V_sub, bsize, k + Param::NUM_COVS);
      Init(V_mask_sub, bsize, k + Param::NUM_COVS);

      Mat<ZZ_p> g, g_mask;
      Vec<ZZ_p> miss, miss_mask;
      g.SetDims(3, m2);
      miss.SetLength(m2);
      g_mask.SetDims(3, m2);
      miss_mask.SetLength(m2);

      Vec<ZZ_p> inner_sx, inner_sxx, inner_sxp;
      Mat<ZZ_p> inner_B;
      Init(inner_sx, m2);
      Init(inner_sxx, m2);
      Init(inner_sxp, m2);
      Init(inner_B, k + Param::NUM_COVS, m2);

      ifstream inner_ifs;
      inner_ifs.open(cache(pid, dataset_idx, "input_geno").c_str(), ios::binary);
      if (pid > 0) {
        mpc.ImportSeed(10, inner_ifs);
      } else {
        for (int p = 1; p <= 2; p++) {
          mpc.ImportSeed(10 + p, inner_ifs);
        }
      }

      // avoid error by re-setting the modulus within each thread
      ZZ base_p = conv<ZZ>(Param::BASE_P.c_str());
      ZZ_p::init(base_p);

      long offset = 0;
      long inner_ind = -1;
      for (int j = 0; j < dataset_idx; j++) {
        offset += n1_vec[j];
        inner_ind += Param::NUM_INDS[j];
      }
      long inner_n1 = n1_vec[dataset_idx];
      for (int cur = 0; cur < inner_n1; cur++) {
        inner_ind++;

        Mat<ZZ_p> g0, g0_mask;
        Vec<ZZ_p> miss0, miss0_mask;

        while (ikeep[inner_ind] != 1) {
          if (pid > 0) {
            mpc.SkipData(inner_ifs, 3, m0); // g
            mpc.SkipData(inner_ifs, m0); // miss

            mpc.SwitchSeed(10);
            mpc.RandMat(g0_mask, 3, m0);
            mpc.RandVec(miss0_mask, m0);
            mpc.RestoreSeed();
          } else {
            for (int p = 1; p <= 2; p++) {
              mpc.SwitchSeed(10 + p);
              mpc.RandMat(g0_mask, 3, m0);
              mpc.RandVec(miss0_mask, m0);
              mpc.RestoreSeed();
            }
          }
          inner_ind++;
        }

        if (pid > 0) {
          mpc.ReadFromFile(g0, inner_ifs, 3, m0); // g
          mpc.ReadFromFile(miss0, inner_ifs, m0); // miss

          mpc.SwitchSeed(10);
          mpc.RandMat(g0_mask, 3, m0);
          mpc.RandVec(miss0_mask, m0);
          mpc.RestoreSeed();
        } else {
          Init(g0, 3, m0);
          Init(g0_mask, 3, m0);
          Init(miss0, m0);
          Init(miss0_mask, m0);
          Vec<ZZ_p> rand_vec;
          Mat<ZZ_p> rand_mat;

          for (int p = 1; p <= 2; p++) {
            mpc.SwitchSeed(10 + p);
            mpc.RandMat(rand_mat, 3, m0);
            mpc.RandVec(rand_vec, m0);
            mpc.RestoreSeed();

            g0_mask += rand_mat;
            miss0_mask += rand_vec;
          }
        }
      
        int inner_ind2 = 0;
        for (int j = 0; j < m0; j++) {
          if (gkeep3[j]) {
            for (int k = 0; k < 3; k++) {
              g[k][inner_ind2] = g0[k][j];
              g_mask[k][inner_ind2] = g0_mask[k][j];
            }
            miss[inner_ind2] = miss0[j];
            miss_mask[inner_ind2] = miss0_mask[j];
            inner_ind2++;
          }
        }

        dosage[cur % bsize] = g[1] + 2 * g[2];
        dosage_mask[cur % bsize] = g_mask[1] + 2 * g_mask[2];

        u_vec[cur % bsize] = u[cur + offset];
        u_mask_vec[cur % bsize] = u_mask[cur + offset];
        p_hat_vec[cur % bsize] = p_hat[cur + offset];
        p_hat_mask_vec[cur % bsize] = p_hat_mask[cur + offset];

        V_sub[cur % bsize] = V[cur + offset];
        V_mask_sub[cur % bsize] = V_mask[cur + offset];

        long new_bsize = bsize;
        if (cur % bsize == bsize - 1  || cur == inner_n1 - 1) {
          if (cur % bsize < bsize - 1) {
            new_bsize = inner_n1 % bsize;
            dosage.SetDims(new_bsize, m2);
            dosage_mask.SetDims(new_bsize, m2);
            u_vec.SetLength(new_bsize);
            u_mask_vec.SetLength(new_bsize);
            p_hat_vec.SetLength(new_bsize);
            p_hat_mask_vec.SetLength(new_bsize);
            V_sub.SetDims(new_bsize, k + Param::NUM_COVS);
            V_mask_sub.SetDims(new_bsize, k + Param::NUM_COVS);
          }

          mpc.ProfilerPopState(false); // file_io/rng

          mpc.BeaverMult(inner_sx, u_vec, u_mask_vec, dosage, dosage_mask);
          mpc.BeaverMult(inner_sxp, p_hat_vec, p_hat_mask_vec, dosage, dosage_mask);

          Mat<ZZ_p> sxx_tmp;
          Init(sxx_tmp, new_bsize, m2);
          mpc.BeaverMultElem(sxx_tmp, dosage, dosage_mask, dosage, dosage_mask);
          for (int b = 0; b < new_bsize; b++) {
            inner_sxx += sxx_tmp[b];
          }
          sxx_tmp.kill();

          mpc.Transpose(V_sub); // (k + NUM_COVS)-by-bsize
          transpose(V_mask_sub, V_mask_sub);

          mpc.BeaverMult(inner_B, V_sub, V_mask_sub, dosage, dosage_mask);

          cout << "\t" << cur + 1 << " / " << inner_n1 << ", "; toc(); tic();

          Init(dosage, bsize, m2);
          Init(dosage_mask, bsize, m2);
          Init(V_sub, bsize, k + Param::NUM_COVS);
          Init(V_mask_sub, bsize, k + Param::NUM_COVS);

          mpc.ProfilerPushState("file_io/rng");
        }
      }
      inner_ifs.close();

      // Add to running sums - each operation mutates shared data, and therefore must be atomic
      #pragma omp critical ( sx_update )
        sx += inner_sx;
      #pragma omp critical ( sxp_update )
        sxp += inner_sxp;
      #pragma omp critical ( sxx_update )
        sxx += inner_sxx;
      #pragma omp critical ( B_update )
        B += inner_B;
    }
    mpc.ProfilerPopState(false); // file_io/rng

    mpc.BeaverReconstruct(sx);
    mpc.BeaverReconstruct(sxp);
    mpc.BeaverReconstruct(sxx);
    mpc.BeaverReconstruct(B);
    sxx *= fp_one;

    fs.open(cache(pid, "gwas_stats").c_str(), ios::out | ios::binary);
    mpc.WriteToFile(sx, fs);
    mpc.WriteToFile(sxx, fs);
    mpc.WriteToFile(sxp, fs);
    mpc.WriteToFile(B, fs);
    fs.close();

    cout << "Wrote results to cache" << endl;
  }

  mpc.ProfilerPopState(true); // data_scan

  if (Param::DEBUG) {
    cout << "sx" << endl;
    mpc.PrintFP(sx, 3);
    cout << "sxp" << endl;
    mpc.PrintFP(sxp, 3);
    cout << "sxx" << endl;
    mpc.PrintFP(sxx, 3);
    cout << "B" << endl;
    mpc.PrintFP(B, 3, 3);
  }

  mpc.Transpose(B); // m2-by-(k + Param::NUM_COVS)

  Vec<ZZ_p> BB;
  mpc.InnerProd(BB, B); // m2
  mpc.FastTrunc(BB);
  if (pid > 0) {
    sxx -= BB;
  }

  ZZ_p sp(0);
  if (pid > 0) {
    for (int i = 0; i < n1; i++) {
      sp += p_hat_mask[i];
      if (pid == 1) {
        sp += p_hat[i];
      }
    }
  }

  ZZ_p spp(0);
  mpc.BeaverInnerProd(spp, p_hat, p_hat_mask);
  mpc.BeaverReconstruct(spp);

  ZZ_p fp_n1_inv = DoubleToFP(1 / ((double) n1), Param::NBIT_K, Param::NBIT_F);
  sx *= fp_n1_inv;
  sp *= fp_n1_inv;

  mpc.FastTrunc(sx);
  mpc.Trunc(sp);
  mpc.Trunc(spp);

  Vec<ZZ_p> sx_mask;
  mpc.BeaverPartition(sx_mask, sx);

  ZZ_p sp_mask;
  mpc.BeaverPartition(sp_mask, sp);

  Vec<ZZ_p> spsx, sx2;
  ZZ_p sp2(0);
  Init(spsx, m2);
  Init(sx2, m2);

  mpc.BeaverMult(spsx, sx, sx_mask, sp, sp_mask);
  mpc.BeaverMult(sp2, sp, sp_mask, sp, sp_mask);
  mpc.BeaverMultElem(sx2, sx, sx_mask, sx, sx_mask);

  mpc.BeaverReconstruct(spsx);
  mpc.BeaverReconstruct(sp2);
  mpc.BeaverReconstruct(sx2);

  spsx *= n1;
  sp2 *= n1;
  sx2 *= n1;

  mpc.FastTrunc(spsx);
  mpc.Trunc(sp2);
  mpc.FastTrunc(sx2);

  Vec<ZZ_p> numer, denom;
  Init(numer, m2);
  Init(denom, m2 + 1);
  if (pid > 0) {
    numer = sxp - spsx;
    for (int i = 0; i < m2; i++) {
      denom[i] = sxx[i] - sx2[i];
    }
    denom[m2] = spp - sp2;
  }

  Vec<ZZ_p> denom1_sqrt_inv;
  if (exists(cache(pid, "denom_inv"))) {
    cout << "denom_inv cache found" << endl;
    ifs.open(cache(pid, "denom_inv").c_str(), ios::binary);
    mpc.ReadFromFile(denom1_sqrt_inv, ifs, denom.length());
    ifs.close();
  } else {
    cout << "Begin denom_inv calculation ... " << endl; tic();
    mpc.ProfilerPushState("sqrt");
    mpc.FPSqrtParallel(tmp_vec, denom1_sqrt_inv, denom);
    mpc.ProfilerPopState(false); // sqrt
    cout << "done. "; toc();

    fs.open(cache(pid, "denom_inv").c_str(), ios::out | ios::binary);
    if (pid > 0) {
      mpc.WriteToFile(denom1_sqrt_inv, fs);
    }
    fs.close();
  }

  denom.kill();
  tmp_vec.kill();

  ZZ_p denom2_sqrt_inv = denom1_sqrt_inv[m2]; // p term
  denom1_sqrt_inv.SetLength(m2); // truncate

  Vec<ZZ_p> z;
  mpc.MultElem(z, numer, denom1_sqrt_inv);
  mpc.FastTrunc(z);

  mpc.MultMat(z, z, denom2_sqrt_inv);
  mpc.FastTrunc(z);

  mpc.ProfilerPopState(false); // assoc_test

  cout << "Association statistics calculated" << endl;
  mpc.RevealSym(z);
  if (pid == 2) {
    Vec<double> z_double;
    cout << "Converting association statistics to double ... "; tic();
    FPToDouble(z_double, z, Param::NBIT_K, Param::NBIT_F);
    cout << "done. "; toc();
    ofs.open(outname("assoc").c_str(), ios::out);
    for (int i = 0; i < z_double.length(); i++) {
      ofs << z_double[i] << endl;
    }
    ofs.close();
    cout << "Result written to " << outname("assoc") << endl;
  }

  mpc.ProfilerPopState(true); // main

  return true;
}

#endif
