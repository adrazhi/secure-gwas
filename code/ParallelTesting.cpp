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
        cout << vect[i];
    }
    cout << endl;
}

int main(int argc, char** argv) {
  if (argc < 5) {
    cout << "Usage: ParallelTesting party_id param_file num_elements num_threads" << endl;
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

  string n_str(argv[3]);
  int n = stoi(n_str);
  cout << "Number of elements in array: " << n << endl;
  
  string num_str(argv[4]);
  int num_threads = stoi(num_str);
  Param::NUM_THREADS = num_threads;
  cout << "Number of threads: " << Param::NUM_THREADS << endl;

  vector< pair<int, int> > pairs;
  pairs.push_back(make_pair(0, 1));
  pairs.push_back(make_pair(0, 2));
  pairs.push_back(make_pair(1, 2));
  MPCEnv mpc;
  if (!mpc.Initialize(pid, pairs)) {
    cout << "MPC environment initialization failed" << endl;
    return 1;
  }

  if (pid == 0) {
    mpc.ReceiveBool(2);
  } else if (pid == 1) {
    // Reconstruct the random mask
    Vec<ZZ_p> a, b;
    mpc.SwitchSeed(2);
    MPCEnv::RandVec(a, n);
    MPCEnv::RandVec(b, n);
    mpc.RestoreSeed();

    // Divide serially
    Vec<ZZ_p> c1;
    mpc.FPDiv(c1, a, b);

    // Divide in parallel
    Vec<ZZ_p> c2;
    mpc.FPDivParallel(c2, a, b);
  } else {
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
    Vec<ZZ_p> a, b;
    Init(a, n);
    Init(b, n);
    for (int i = 0; i < n; i++) {
      DoubleToFP(a[i], a_base[i], Param::NBIT_K, Param::NBIT_F);
      DoubleToFP(b[i], b_base[i], Param::NBIT_K, Param::NBIT_F);
    }
    
    // Now generate the random mask and mask out the data
    Vec<ZZ_p> ra, rb;
    mpc.SwitchSeed(1);
    MPCEnv::RandVec(ra, n);
    MPCEnv::RandVec(rb, n);
    mpc.RestoreSeed();
    a -= ra;
    b -= rb;

    // Party 2 prints profiling, i.e. runtime, data
    struct timeval start, end;
    double runtime;
    
    // Divide serially
    Vec<ZZ_p> c1;
    gettimeofday(&start, NULL); 
    ios_base::sync_with_stdio(false);
    mpc.FPDiv(c1, a, b);
    gettimeofday(&end, NULL);

    Vec<double> c1_base;
    FPToDouble(c1_base, c1, Param::NBIT_K, Param::NBIT_F);
    print_vec("Division (serial)", c1_base, 5);
    runtime = (end.tv_sec - start.tv_sec) * 1e6;
    runtime = (runtime + (end.tv_usec - start.tv_usec)) * 1e-6;
    cout << "Runtime (serial): " << fixed << runtime << setprecision(6); 
    cout << " sec" << endl;

    // Divide in parallel
    Vec<ZZ_p> c2;
    gettimeofday(&start, NULL); 
    ios_base::sync_with_stdio(false);
    mpc.FPDivParallel(c2, a, b);
    gettimeofday(&end, NULL);

    Vec<double> c2_base;
    FPToDouble(c2_base, c2, Param::NBIT_K, Param::NBIT_F);
    print_vec("Division (parallel)", c2_base, 5);
    runtime = (end.tv_sec - start.tv_sec) * 1e6;
    runtime = (runtime + (end.tv_usec - start.tv_usec)) * 1e-6;
    cout << "Runtime (parallel): " << fixed << runtime << setprecision(6); 
    cout << " sec" << endl;

    // All done, so send signal to party 0
    mpc.SendBool(true, 0);
  }

  mpc.CleanUp();

  cout << "Protocol successfully completed" << endl;
  return 0;
}
