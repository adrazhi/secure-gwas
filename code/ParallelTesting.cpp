#include "connect.h"
#include "mpc.h"
#include "protocol.h"
#include "util.h"
#include "NTL/ZZ_p.h"
#include "gwasiter.h"

#include <bits/stdc++.h>
#include <sys/time.h>
#include <cstdlib>
#include <fstream>
#include <map>
#include <iostream>
#include <string>
#include <sstream>
#include <random>
#include <vector>
#include <chrono>

using namespace NTL;
using namespace std;

using msec = chrono::milliseconds;
using get_time = chrono::steady_clock;

void print_vec(const char* name, vector<double> &vect, int num) {
  cout << name << ": ";
  for (int i = 0; i < num; i++) {
    cout << vect[i] << "  ";
  }
  cout << endl;
}

int main(int argc, char** argv) {
  if (argc < 4) {
    cout << "Usage: ParallelTesting party_id param_file num_elements" << endl;
    return 1;
  }

  string pid_str(argv[1]);
  int pid;
  if (!Param::Convert(pid_str, pid, "party_id") || pid < 0 || pid > 2) {
    cout << "Error: party_id should be 0, 1, or 2" << endl;
    return 1;
  }

  if (!Param::ParseFile(argv[2])) {
    cout << "Could not finish parsing parameter file" << endl;
    return 1;
  }

  vector<int> num_threads{ 1, 2, 4, 8, 16, 32};
  Param::NUM_THREADS = 32; // need this so that mpc.Initialize will create enough channels/prgs for max number of threads

  string n_str(argv[3]);
  long n = stoi(n_str);
  cout << "Number of elements in array: " << n << endl;

  int n_trials = 5;
  if (argc == 5) {
    string n_tr(argv[4]);
    n_trials = stoi(n_tr);
  }

  vector< pair<int, int> > pairs;
  pairs.push_back(make_pair(0, 1));
  pairs.push_back(make_pair(0, 2));
  pairs.push_back(make_pair(1, 2));
  MPCEnv mpc;
  if (!mpc.Initialize(pid, pairs)) {
    cout << "MPC environment initialization failed" << endl;
    return 1;
  }
  // mpc.SetDebug(true);

  // Initialize vectors to hold inputs and outputs of parallel subroutines
  Vec<ZZ_p> a, b, c1, c2;
  Init(a, n);
  Init(b, n);
  Init(c1, n);
  Init(c2, n);

  Mat<ZZ_p> Y, Q;
  Init(Y, 15, n);

  if (pid == 1) {
    // Reconstruct the random mask
    mpc.SwitchSeed(2);
    mpc.RandVec(a, n);
    mpc.RandVec(b, n);
    mpc.RandMat(Y, 15, n);
    mpc.RestoreSeed();
  } else if (pid == 2) {
    // Generate vectors of random doubles to simulate data
    std::uniform_real_distribution<double> unif(1, 10);
    std::default_random_engine re;
    auto rand_dbl = std::bind(unif, re);
    vector<double> a_base, b_base;
    for (int i = 0; i < n; i++) {
      a_base.push_back(rand_dbl());
      b_base.push_back(rand_dbl());
    }
    print_vec("Vector 1", a_base, 5);
    print_vec("Vector 2", b_base, 5);

    // Convert the data from double to FP
    cout << "Converting double to FP ... ";
    for (int i = 0; i < n; i++) {
      DoubleToFP(a[i], a_base[i], Param::NBIT_K, Param::NBIT_F);
      DoubleToFP(b[i], b_base[i], Param::NBIT_K, Param::NBIT_F);
    }
    cout << "done" << endl;

    // Generate random matrix
    mpc.RandMat(Y, 15, n);
    
    // Now generate the random mask and mask out the data
    cout << "Masking data ... ";
    Vec<ZZ_p> ra, rb;
    Mat<ZZ_p> rY;
    mpc.SwitchSeed(1);
    mpc.RandVec(ra, n);
    mpc.RandVec(rb, n);
    mpc.RandMat(rY, 15, n);
    mpc.RestoreSeed();
    a -= ra;
    b -= rb;
    Y -= rY;
    cout << "done" << endl;
  }

  // output flags
  bool fpdiv = true;
  bool fpsqrt = true;
  bool ispos = true;
  bool orth = true;
  bool print_output = false;

  for (int i = 0; i < num_threads.size(); i++) {
    Param::NUM_THREADS = num_threads[i];
    cout << "-----------------" << endl;
    cout << "Number of Threads: " << Param::NUM_THREADS << endl;

    // Profile FPDivParallel
    if (fpdiv) {
      auto start = get_time::now();
      for (int i = 0; i < n_trials; i++) {
        mpc.FPDivParallel(c1, a, b);
      }
      auto end = get_time::now();

      // Print results
      double runtime = chrono::duration_cast<msec>(end - start).count() / (n_trials * 1000.0);
      cout << "Avg Division Runtime: " << runtime << " sec" << endl;
      if (print_output) {
        cout << "Division Result: ";
        mpc.PrintFP(c1, 5);
      }
    }

    // Profile FPSqrtParallel
    if (fpsqrt) {
      auto start = get_time::now();
      for (int i = 0; i < n_trials; i++) {
        mpc.FPSqrtParallel(c1, c2, a);
      }
      auto end = get_time::now();

      // Print results
      double runtime = chrono::duration_cast<msec>(end - start).count() / (n_trials * 1000.0);
      cout << "Avg Sqrt Runtime: " << runtime << " sec" << endl;
      if (print_output) {
        cout << "Sqrt Result 1: ";
        mpc.PrintFP(c1, 5);
        cout << "Sqrt Result 2: ";
        mpc.PrintFP(c2, 5);
      }
    }

    // Profile IsPositiveParallel
    if (ispos) {
      Vec<ZZ_p> input = a - b;
      auto start = get_time::now();
      for (int i = 0; i < n_trials; i++) {
        mpc.IsPositiveParallel(c1, input);
      }
      auto end = get_time::now();

      // Print results
      double runtime = chrono::duration_cast<msec>(end - start).count() / (n_trials * 1000.0);
      cout << "Avg IsPositive Runtime: " << runtime << " sec" << endl;
      if (print_output) {
        cout << "IsPositive Result: ";
        mpc.PrintFP(c1, 5);
      }
    }

    // Profile OrthonormalBasis
    if (orth) {
      auto start = get_time::now();
      for (int i = 0; i < n_trials; i++) {
        mpc.OrthonormalBasis(Q, Y);
      }
      auto end = get_time::now();

      // Print results
      double runtime = chrono::duration_cast<msec>(end - start).count() / (n_trials * 1000.0);
      cout << "Avg OrthonormalBasis Runtime: " << runtime << " sec" << endl;
    }
  }

  mpc.CleanUp();

  cout << "Protocol successfully completed" << endl;
  return 0;
}
