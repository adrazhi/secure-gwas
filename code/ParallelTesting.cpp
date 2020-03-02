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
#include <sstream>

using namespace NTL;
using namespace std;

int main(int argc, char** argv) {
  if (!Param::ParseFile(argv[1])) {
    cout << "Could not finish parsing parameter file" << endl;
    return 1;
  }

  int n = 100000;
  Param::NUM_THREADS = 10;

  Vec<ZZ_p> a, b, c1, c2;
  MPCEnv::RandVec(a, n);
  MPCEnv::RandVec(b, n);
  struct timeval start, end;

  gettimeofday(&start, NULL); 
  ios_base::sync_with_stdio(false); 
  MPCEnv::FPDiv(c1, a, b);
  gettimeofday(&end, NULL); 

  double runtime;
  runtime = (end.tv_sec - start.tv_sec) * 1e6;
  runtime = (runtime + (end.tv_usec - start.tv_usec)) * 1e-6;
  cout << "Runtime (serial): " << fixed << time_taken << setprecision(6); 
  cout << " sec" << endl;

  gettimeofday(&start, NULL); 
  ios_base::sync_with_stdio(false); 
  MPCEnv::FPDivParallel(c2, a, b);
  gettimeofday(&end, NULL); 

  runtime = (end.tv_sec - start.tv_sec) * 1e6;
  runtime = (runtime + (end.tv_usec - start.tv_usec)) * 1e-6;
  cout << "Runtime (parallel): " << fixed << time_taken << setprecision(6); 
  cout << " sec" << endl; 
}