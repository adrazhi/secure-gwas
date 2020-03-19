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

using namespace NTL;
using namespace std;

int main(int argc, char** argv) {
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

  Vec<ZZ_p> a, b, c1, c2;
  mpc.SwitchSeed(0);
  mpc.RandVec(a, n);
  cout << "Vector 1: " << a[0] << endl;
  mpc.RandVec(b, n);
  cout << "Vector 2: " << b[0] << endl;
  mpc.RestoreSeed();
  struct timeval start, end;
  double runtime;

  gettimeofday(&start, NULL); 
  ios_base::sync_with_stdio(false);
  mpc.ProfilerPushState("div"); 
  mpc.FPDiv(c1, a, b);
  mpc.ProfilerPopState(false); // div
  gettimeofday(&end, NULL); 

  cout << "Division 1: " << c1[0] << endl;

  runtime = (end.tv_sec - start.tv_sec) * 1e6;
  runtime = (runtime + (end.tv_usec - start.tv_usec)) * 1e-6;
  cout << "Runtime (serial): " << fixed << runtime << setprecision(6); 
  cout << " sec" << endl;

  gettimeofday(&start, NULL); 
  ios_base::sync_with_stdio(false); 
  mpc.FPDivParallel(c2, a, b);
  gettimeofday(&end, NULL); 

  cout << "Division 2: " << c2[0] << endl;

  runtime = (end.tv_sec - start.tv_sec) * 1e6;
  runtime = (runtime + (end.tv_usec - start.tv_usec)) * 1e-6;
  cout << "Runtime (parallel): " << fixed << runtime << setprecision(6); 
  cout << " sec" << endl; 

  mpc.CleanUp();
}