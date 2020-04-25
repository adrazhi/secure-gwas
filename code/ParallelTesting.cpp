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

using namespace NTL;
using namespace std;

void print_vec(const char* name, vector<double> &vect, int num) {
  cout << name << ": ";
  for (int i = 0; i < num; i++) {
    cout << vect[i] << "  ";
  }
  cout << endl;
}

void print_ntl_vec(const char* name, Vec<double> &vect, int num) {
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

  vector<int> num_threads{ 1, 2, 4, 8, 16 };
  Param::NUM_THREADS = 16; // need this so that mpc.Initialize will create enough channels/prgs for up to 16 threads

  string n_str(argv[3]);
  long n = stoi(n_str);
  cout << "Number of elements in array: " << n << endl; 

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
  Vec<double> c1_base;
  Vec<double> c2_base;

  if (pid == 1) {
    // Reconstruct the random mask
    mpc.SwitchSeed(2);
    mpc.RandVec(a, n);
    mpc.RandVec(b, n);
    mpc.RestoreSeed();
  } else if (pid == 2) {
    // Generate vectors of random doubles to simulate data
    std::uniform_real_distribution<double> unif(-10, 10);
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
    
    // Now generate the random mask and mask out the data
    cout << "Masking data ... ";
    Vec<ZZ_p> ra, rb;
    mpc.SwitchSeed(1);
    mpc.RandVec(ra, n);
    mpc.RandVec(rb, n);
    mpc.RestoreSeed();
    a -= ra;
    b -= rb;
    cout << "done" << endl;
  }

  struct timeval start, end;
  double runtime;

  // output flags
  bool fpdiv = true;
  bool fpsqrt = true;
  bool ispos = true;
  bool print_output = true;

  for (int i = 0; i < num_threads.size(); i++) {
    Param::NUM_THREADS = num_threads[i];
    cout << "-----------------" << endl;
    cout << "Number of Threads: " << Param::NUM_THREADS << endl;

    // Profile FPDivParallel
    if (fpdiv) {
      gettimeofday(&start, NULL); 
      ios_base::sync_with_stdio(false);
      mpc.FPDivParallel(c1, a, b);
      gettimeofday(&end, NULL);
      mpc.RevealSym(c1);

      // Print results
      runtime = (end.tv_sec - start.tv_sec) * 1e6;
      runtime = (runtime + (end.tv_usec - start.tv_usec)) * 1e-6;
      FPToDouble(c1_base, c1, Param::NBIT_K, Param::NBIT_F);
      if (pid == 2) {
        cout << "Division Runtime: " << fixed << runtime << setprecision(6); 
        cout << " sec" << endl;
        if (print_output) print_ntl_vec("Division Result", c1_base, 5);
        cout << endl;
      }
    }

    // Profile FPSqrtParallel
    if (fpsqrt) {
      gettimeofday(&start, NULL); 
      ios_base::sync_with_stdio(false);
      mpc.FPSqrtParallel(c1, c2, a);
      gettimeofday(&end, NULL);
      mpc.RevealSym(c1);
      mpc.RevealSym(c2);

      // Print results
      runtime = (end.tv_sec - start.tv_sec) * 1e6;
      runtime = (runtime + (end.tv_usec - start.tv_usec)) * 1e-6;
      FPToDouble(c1_base, c1, Param::NBIT_K, Param::NBIT_F);
      FPToDouble(c2_base, c2, Param::NBIT_K, Param::NBIT_F);
      if (pid == 2) {
        cout << "Square Root Runtime: " << fixed << runtime << setprecision(6); 
        cout << " sec" << endl;
        if (print_output) {
          print_ntl_vec("Square Root Result 1", c1_base, 5);
          print_ntl_vec("Square Root Result 2", c2_base, 5);
        }
      }
    }

    // Profile IsPositiveParallel
    if (fpdiv) {
      gettimeofday(&start, NULL); 
      ios_base::sync_with_stdio(false);
      mpc.IsPositiveParallel(c1, a);
      gettimeofday(&end, NULL);
      mpc.RevealSym(c1);

      // Print results
      runtime = (end.tv_sec - start.tv_sec) * 1e6;
      runtime = (runtime + (end.tv_usec - start.tv_usec)) * 1e-6;
      if (pid == 2) {
        cout << "IsPositive Runtime: " << fixed << runtime << setprecision(6); 
        cout << " sec" << endl;
      }
      if (print_output) {
        mpc.Print(c1, 5);
      }
    }
  }

  mpc.CleanUp();

  cout << "Protocol successfully completed" << endl;
  return 0;
}
