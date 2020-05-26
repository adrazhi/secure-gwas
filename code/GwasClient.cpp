#include "connect.h"
#include "mpc.h"
#include "protocol.h"
#include "util.h"
#include "NTL/ZZ_p.h"
#include "gwasiter.h"

#include <cstdlib>
#include <fstream>
#include <map>
#include <string>
#include <iostream>
#include <sstream>
#include <bits/stdc++.h>
#include <sys/time.h>

using namespace NTL;
using namespace std;

// simple helper method to print a vector
template <class T>
void print_vec(const char* name, vector<T> &vect, int num) {
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

  // reset the round number, since all datasets will be joined during GWAS execution
  Param::CUR_ROUND = 0;

  if (!Param::ParseFile(argv[2])) {
    cout << "Could not finish parsing parameter file" << endl;
    return 1;
  }

  // pre-process the param class before running GWAS
  // in particular, we want to create a list of cache file prefixes for every chunk that was created during Data Sharing,
  // and we want to update the NUM_INDS vector correspondingly so that NUM_INDS[i] = size(chunk[i])
  int n = Param::NUM_INDS.size();
  int total_chunks = 0;
  // iterate over all datasets
  for (int i = 0; i < n; i++) {
    int num_chunks = Param::NUM_CHUNKS[i];
    total_chunks += num_chunks;

    // replicate the chunk division that was used to split datasets during Data Sharing
    long total = Param::NUM_INDS[i];
    long chunk_size = ceil(total / ((double) num_chunks));
    num_chunks = ceil(total / ((double) chunk_size)); 
    long remainder = total - ((num_chunks - 1) * chunk_size);

    // append each chunk size to the end of the NUM_INDS vector
    for (int j = 0; j < num_chunks - 1; j++) {
      Param::NUM_INDS.push_back(chunk_size);
    }
    Param::NUM_INDS.push_back(remainder);

    // append each chunk file prefix to the end of the CACHE_FILE_PREFIX vector
    for (int j = 0; j < num_chunks; j++) {
      string pref = Param::CACHE_FILE_PREFIX[i] + "_" + to_string(j);
      Param::CACHE_FILE_PREFIX.push_back(pref);
    }
  }

  // now get rid of old values that correspond to dataset information instead of chunk information
  for (int i = n; i < n + total_chunks; i++) {
    Param::NUM_INDS[i-n] = Param::NUM_INDS[i];
    Param::CACHE_FILE_PREFIX[i-n] = Param::CACHE_FILE_PREFIX[i];
  }
  Param::NUM_INDS.resize(total_chunks);
  Param::CACHE_FILE_PREFIX.resize(total_chunks);

  print_vec("NUM_INDS", Param::NUM_INDS, Param::NUM_INDS.size());
  print_vec("CACHE_FILE_PREFIX", Param::CACHE_FILE_PREFIX, Param::CACHE_FILE_PREFIX.size());

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
