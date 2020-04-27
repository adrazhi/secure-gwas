#include "connect.h"
#include "mpc.h"
#include "protocol.h"
#include "util.h"
#include "NTL/ZZ_p.h"
#include "gwasiter.h"
#include <NTL/BasicThreadPool.h>

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

int main(int argc, char** argv) {
  if (argc < 3) {
    cout << "Usage: Benchmarking party_id param_file" << endl;
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

  Param::NUM_THREADS = 20;
  cout << "Num Threads: " << Param::NUM_THREADS << endl;

  // SetNumThreads(5);
  // cout << AvailableThreads() << " threads created for NTL" << endl;

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

  // Initialize matrices to hold inputs and outputs
  Mat<ZZ_p> Y1, Y2, Y3, Q;
  Init(Y1, 15, 1000);
  Init(Y2, 15, 10000);
  Init(Y3, 15, 100000);

  Mat<ZZ_p> C, A, B;
  Init(A, 15, 1);
  Init(B, 1, 100000);

  if (pid == 1) {
    // Reconstruct the random mask
    mpc.SwitchSeed(2);
    mpc.RandMat(Y1, 15, 1000);
    mpc.RandMat(Y2, 15, 10000);
    mpc.RandMat(Y3, 15, 100000);

    mpc.RandMat(A, 15, 1);
    mpc.RandMat(B, 1, 100000);
    mpc.RestoreSeed();
  } else if (pid == 2) {
    // Generate data
    mpc.RandMat(Y1, 15, 1000);
    mpc.RandMat(Y2, 15, 10000);
    mpc.RandMat(Y3, 15, 100000);

    mpc.RandMat(A, 15, 1);
    mpc.RandMat(B, 1, 100000);
    
    // Mask out data
    cout << "Masking data ... ";
    Mat<ZZ_p> r1, r2, r3;
    mpc.SwitchSeed(1);
    mpc.RandMat(r1, 15, 1000);
    mpc.RandMat(r2, 15, 10000);
    mpc.RandMat(r3, 15, 100000);
    Mat<ZZ_p> r4, r5;
    mpc.RandMat(r4, 15, 1);
    mpc.RandMat(r5, 1, 100000);
    mpc.RestoreSeed();
    Y1 -= r1;
    Y2 -= r2;
    Y3 -= r3;
    A -= r4;
    B -= r5;
    cout << "done" << endl;
  }

  tic(); mpc.MultMat(C, A, B); toc();
  tic(); mpc.Trunc(C); toc();
  cout << "----" << endl;
  tic(); mpc.FastMultMat(C, A, B); toc();
  tic(); mpc.FastTrunc(C); toc();
  cout << "----" << endl;
  transpose(B);
  tic(); mpc.MultMat(C, Y3, B); toc();
  tic(); mpc.Trunc(C); toc();
  cout << "----" << endl;
  tic(); mpc.FastMultMat(C, Y3, B); toc();
  tic(); mpc.FastTrunc(C); toc();
  cout << "----" << endl;
  transpose(B);
  tic(); mpc.FastMultMat2(C, A, B); toc();
  tic(); mpc.FastTrunc(C); toc();

  // Param::NUM_THREADS = 1;
  // tic(); mpc.OrthonormalBasis(Q, Y3); toc();
  // cout << "-----------" << endl;
  // Param::NUM_THREADS = 20;
  // tic(); mpc.OrthonormalBasis(Q, Y3); toc();

  mpc.CleanUp();

  cout << "Protocol successfully completed" << endl;
  return 0;
}
