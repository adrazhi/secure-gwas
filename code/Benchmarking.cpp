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
  if (!Param::Convert(pid_str, pid, "party_id") || pid < 0 || pid > 3) {
    cout << "Error: party_id should be 0, 1, 2, or 3" << endl;
    return 1;
  }

  if (!Param::ParseFile(argv[2])) {
    cout << "Could not finish parsing parameter file" << endl;
    return 1;
  }

  Param::NUM_THREADS = 20;
  cout << "Num Threads: " << Param::NUM_THREADS << endl;

  vector< pair<int, int> > pairs;
  pairs.push_back(make_pair(0, 1));
  pairs.push_back(make_pair(0, 2));
  pairs.push_back(make_pair(1, 2));
  // pairs.push_back(make_pair(1, 3));
  // pairs.push_back(make_pair(2, 3));
  MPCEnv mpc;
  if (!mpc.Initialize(pid, pairs)) {
    cout << "MPC environment initialization failed" << endl;
    return 1;
  }
  // mpc.SetDebug(true);

  // OrthonormalBasis Benchmarking

  // // Initialize matrices to hold inputs and outputs
  Mat<ZZ_p> Y1, Y2, Y3, Q;
  Init(Y1, 15, 2000);
  Init(Y2, 15, 20000);
  Init(Y3, 15, 200000);

  // Mat<ZZ_p> C, A, B, B2;
  // Init(A, 15, 1);
  // Init(B, 1, 100000);
  // Init(B2, 100000, 1);


  if (pid == 1) {
    // Reconstruct the random mask
    mpc.SwitchSeed(2);
    mpc.RandMat(Y1, 15, 1000);
    mpc.RandMat(Y2, 15, 10000);
    mpc.RandMat(Y3, 15, 100000);

  //   mpc.RandMat(A, 15, 1);
  //   mpc.RandMat(B, 1, 100000);
  //   mpc.RandMat(B2, 100000, 1);
    mpc.RestoreSeed();
  } else if (pid == 2) {
    // Generate data
    mpc.RandMat(Y1, 15, 1000);
    mpc.RandMat(Y2, 15, 10000);
    mpc.RandMat(Y3, 15, 100000);

  //   mpc.RandMat(A, 15, 1);
  //   mpc.RandMat(B, 1, 100000);
  //   mpc.RandMat(B2, 100000, 1);
    
  //   // Mask out data
    cout << "Masking data ... ";
    Mat<ZZ_p> r1, r2, r3;
    mpc.SwitchSeed(1);
    mpc.RandMat(r1, 15, 1000);
    mpc.RandMat(r2, 15, 10000);
    mpc.RandMat(r3, 15, 100000);
  //   Mat<ZZ_p> r4, r5, r6;
  //   mpc.RandMat(r4, 15, 1);
  //   mpc.RandMat(r5, 1, 100000);
  //   mpc.RandMat(r6, 100000, 1);
    mpc.RestoreSeed();
    Y1 -= r1;
    Y2 -= r2;
    Y3 -= r3;
  //   A -= r4;
  //   B -= r5;
  //   B2 -= r6;
    cout << "done" << endl;
  }

  for (int i = 0; i < 3; i++) { cout << "Trunc (15 x 2,000) ... "; tic(); mpc.Trunc(Y1); toc(); }
  cout << "---" << endl;
  for (int i = 0; i < 3; i++) { cout << "Trunc (15 x 20,000) ... "; tic(); mpc.Trunc(Y2); toc(); }
  cout << "---" << endl;
  for (int i = 0; i < 3; i++) { cout << "Trunc (15 x 200,000) ... "; tic(); mpc.Trunc(Y3); toc(); }
  cout << "-----------" << endl;
  for (int i = 0; i < 3; i++) { cout << "FastTrunc (15 x 2,000) ... "; tic(); mpc.FastTrunc(Y1); toc(); }
  cout << "---" << endl;
  for (int i = 0; i < 3; i++) { cout << "FastTrunc (15 x 20,000) ... "; tic(); mpc.FastTrunc(Y2); toc(); }
  cout << "---" << endl;
  for (int i = 0; i < 3; i++) { cout << "FastTrunc (15 x 200,000) ... "; tic(); mpc.FastTrunc(Y3); toc(); }
  // // cout << "Serial: Mult (15 by 1) x (1 x 100k) ... "; tic(); mpc.MultMat(C, A, B); toc();
  // // cout << "Serial: Trunc 15 by 100k ... "; tic(); mpc.Trunc(C); toc();
  // // cout << "----" << endl;
  // // cout << "Parallel Mult (15 by 1) x (1 x 100k) ... ";tic(); mpc.FastMultMat(C, A, B); toc();
  // // cout << "Parallel Trunc 15 by 100k ... "; tic(); mpc.FastTrunc(C); toc();
  // // cout << "----" << endl;

  // // cout << "Serial: Mult (15 by 100k) x (100k x 1) ... "; tic(); mpc.MultMat(C, Y3, B2); toc();
  // // cout << "Serial: Trunc 15 by 1 ... "; tic(); mpc.Trunc(C); toc();
  // // cout << "----" << endl;
  // // cout << "Parallel: Mult (15 by 100k) x (100k x 1) ... "; tic(); mpc.FastMultMat(C, Y3, B2); toc();
  // // cout << "Parallel: Trunc 15 by 1 ... "; tic(); mpc.FastTrunc(C); toc();

  // // Param::NUM_THREADS = 1;
  // // tic(); mpc.OrthonormalBasis(Q, Y1); toc();
  // // tic(); mpc.OrthonormalBasis(Q, Y2); toc();
  // // tic(); mpc.OrthonormalBasis(Q, Y3); toc();
  // // cout << "-----------" << endl;
  // // Param::NUM_THREADS = 20;
  // // tic(); mpc.OrthonormalBasis(Q, Y1); toc();
  // // tic(); mpc.OrthonormalBasis(Q, Y2); toc();
  // // tic(); mpc.OrthonormalBasis(Q, Y3); toc();

  // Param::NUM_THREADS = 1;
  // for (int i = 0; i < 10; i++) {
  //   tic(); mpc.OrthonormalBasis(Q, Y1); toc();
  // }
  // Param::NUM_THREADS = 20;
  // for (int i = 0; i < 10; i++) {
  //   tic(); mpc.OrthonormalBasis(Q, Y1); toc();
  // }

  // // Data Transfer Benchmarking
  // int n = 500000;

  // // Between 0 and 1 and between 2 and 3
  // if (pid == 0) {
  //   Vec<ZZ_p> vec;
  //   mpc.RandVec(vec, n);
  //   cout << "Sending to 1" << endl;
  //   for (int i = 0; i < 5; i++) {
  //     tic(); mpc.SendVec(vec, 1); toc();
  //   }
  //   cout << "Receiving from 1" << endl;
  //   for (int i = 0; i < 5; i++) {
  //     tic(); mpc.ReceiveVec(vec, 1, n); toc();
  //   }
  // } else if (pid == 1) {
  //   Vec<ZZ_p> vec;
  //   mpc.RandVec(vec, n);
  //   cout << "Receiving from 0" << endl;
  //   for (int i = 0; i < 5; i++) {
  //     tic(); mpc.ReceiveVec(vec, 0, n); toc();
  //   }
  //   cout << "Sending to 0" << endl;
  //   for (int i = 0; i < 5; i++) {
  //     tic(); mpc.SendVec(vec, 0); toc();
  //   }
  // } else if (pid == 2) {
  //   Vec<ZZ_p> vec;
  //   mpc.RandVec(vec, n);
  //   cout << "Sending to 3" << endl;
  //   for (int i = 0; i < 5; i++) {
  //     tic(); mpc.SendVec(vec, 3); toc();
  //   }
  //   cout << "Receiving from 3" << endl;
  //   for (int i = 0; i < 5; i++) {
  //     tic(); mpc.ReceiveVec(vec, 3, n); toc();
  //   }
  // } else if (pid == 3) {
  //   Vec<ZZ_p> vec;
  //   mpc.RandVec(vec, n);
  //   cout << "Receiving from 2" << endl;
  //   for (int i = 0; i < 5; i++) {
  //     tic(); mpc.ReceiveVec(vec, 2, n); toc();
  //   }
  //   cout << "Sending to 2" << endl;
  //   for (int i = 0; i < 5; i++) {
  //     tic(); mpc.SendVec(vec, 2); toc();
  //   }
  // }

  // // Between 0 and 2 and between 1 and 3
  // if (pid == 0) {
  //   Vec<ZZ_p> vec;
  //   mpc.RandVec(vec, n);
  //   cout << "Sending to 2" << endl;
  //   for (int i = 0; i < 5; i++) {
  //     tic(); mpc.SendVec(vec, 2); toc();
  //   }
  //   cout << "Receiving from 2" << endl;
  //   for (int i = 0; i < 5; i++) {
  //     tic(); mpc.ReceiveVec(vec, 2, n); toc();
  //   }
  // } else if (pid == 2) {
  //   Vec<ZZ_p> vec;
  //   mpc.RandVec(vec, n);
  //   cout << "Receiving from 0" << endl;
  //   for (int i = 0; i < 5; i++) {
  //     tic(); mpc.ReceiveVec(vec, 0, n); toc();
  //   }
  //   cout << "Sending to 0" << endl;
  //   for (int i = 0; i < 5; i++) {
  //     tic(); mpc.SendVec(vec, 0); toc();
  //   }
  // } else if (pid == 1) {
  //   Vec<ZZ_p> vec;
  //   mpc.RandVec(vec, n);
  //   cout << "Sending to 3" << endl;
  //   for (int i = 0; i < 5; i++) {
  //     tic(); mpc.SendVec(vec, 3); toc();
  //   }
  //   cout << "Receiving from 3" << endl;
  //   for (int i = 0; i < 5; i++) {
  //     tic(); mpc.ReceiveVec(vec, 3, n); toc();
  //   }
  // } else if (pid == 3) {
  //   Vec<ZZ_p> vec;
  //   mpc.RandVec(vec, n);
  //   cout << "Receiving from 1" << endl;
  //   for (int i = 0; i < 5; i++) {
  //     tic(); mpc.ReceiveVec(vec, 1, n); toc();
  //   }
  //   cout << "Sending to 1" << endl;
  //   for (int i = 0; i < 5; i++) {
  //     tic(); mpc.SendVec(vec, 1); toc();
  //   }
  // }

  // // Between 1 and 2
  // if (pid == 1) {
  //   Vec<ZZ_p> vec;
  //   mpc.RandVec(vec, n);
  //   cout << "Sending to 2" << endl;
  //   for (int i = 0; i < 5; i++) {
  //     tic(); mpc.SendVec(vec, 2); toc();
  //   }
  //   cout << "Receiving from 2" << endl;
  //   for (int i = 0; i < 5; i++) {
  //     tic(); mpc.ReceiveVec(vec, 2, n); toc();
  //   }
  // } else if (pid == 2) {
  //   Vec<ZZ_p> vec;
  //   mpc.RandVec(vec, n);
  //   cout << "Receiving from 1" << endl;
  //   for (int i = 0; i < 5; i++) {
  //     tic(); mpc.ReceiveVec(vec, 1, n); toc();
  //   }
  //   cout << "Sending to 1" << endl;
  //   for (int i = 0; i < 5; i++) {
  //     tic(); mpc.SendVec(vec, 1); toc();
  //   }
  // }

  mpc.CleanUp();

  cout << "Protocol successfully completed" << endl;
  return 0;
}
