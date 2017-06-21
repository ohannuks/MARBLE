#include <cassert>
#include "cstdlib"
#include <iostream>
#include <chrono> 
#include <ctime>
#include <cstdio>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>
#include "helperFunctions.h"
#include <cmath>
// feenableexcept() 
#if defined HAVE_FENV_H
# include <fenv.h>
#endif
// signal()
#include <csignal>
// getpid()
#include <sys/types.h>
#include <unistd.h>
// MPI
#if defined HAVE_MPI
//#include <mpi.h>
#endif
// ULONG_MAX
#include <climits>
#include <fstream>
#include <random>
#include <cstdlib>
#include <cstdio>
#include <string>
//#include "mpi.h"
#include "omp.h"
#include "sampling.h"


using namespace std;
const int starMaxIterations = 1000000;

const double pi = M_PI;

bool sampleDistribution( double consts[4], mt19937 & mt, const double M, const double a, const double rMinLimit ) {
  // Sample dark matter distributions 
  while( true ) {
    uniform_real_distribution<double> distE (0.9,1.0);
    uniform_real_distribution<double> distLz(-60.0*M,60.0*M);
    uniform_real_distribution<double> distQ(0,(60.*M)*(60.*M));
    double rMin = 0, rMax = 0, thMin = 0, thMax = 0;
    const double Eps = distE(mt), Lz = distLz(mt), Q = distQ(mt);
    const double K   = Q + pow(Eps*a-Lz,2);
    consts[0] = Eps;
    consts[1] = Lz;
    consts[2] = Q;
    consts[3] = K;
    if( helper::kerr_r_limits(Eps, Lz, Q, a, M, rMin, rMax) && helper::kerr_th_limits(Eps, Lz, Q, a, M, thMin, thMax) && rMin < rMinLimit ) { 
      //cout << "Inputting distribution constants with Eps, Lz, Q: " << Eps << " " << Lz << " " << Q << endl;
      break;
    }
  }
  assert( abs(consts[1]) < 52.*M);
  assert( abs(consts[2]) < (52.*M)*(52.*M) ); // If this is true, we might be missing some of the possible orbits
  return true;
}

bool save_results(const vector<double> & EpsVector, 
                  const vector<double> & LzVector,
                  const vector<double> & QVector, 
                  const vector<pair<double,int>> & Sr, 
                  const vector<pair<double,int>> & Sth,
                  const vector<pair<double,int>> & Lz,
                  const int nParticles,
                  const double rMinLimit,
                  const double a,
                  const double M,
                  const char * filename
                 ) {
  assert( EpsVector.size() == nParticles );
  FILE * f = fopen(filename, "wb");

  // Write attributes
  fwrite( (const void*)&nParticles, sizeof(nParticles),    1, f);
  fwrite( (const void*)&rMinLimit,  sizeof(rMinLimit), 1, f);
  fwrite( (const void*)&a,  sizeof(a), 1, f);
  fwrite( (const void*)&M,  sizeof(M), 1, f);
  fwrite( (const void*)EpsVector.data(), sizeof(double), EpsVector.size(), f );
  fwrite( (const void*)LzVector.data(), sizeof(double), LzVector.size(), f );
  fwrite( (const void*)QVector.data(), sizeof(double), QVector.size(), f );
  vector<double> Sr_1 (Sr.size() );
  vector<int>    Sr_2 (Sr.size() );
  vector<double> Sth_1(Sth.size());
  vector<int>    Sth_2(Sth.size());
  vector<double> Lz_1(Lz.size());
  vector<int>    Lz_2(Lz.size());
  for( int i =0; i < Sr.size(); ++i ) {
    Sr_1[i] = Sr[i].first;
    Sr_2[i] = Sr[i].second;
    Sth_1[i]= Sth[i].first;
    Sth_2[i]= Sth[i].second;
    Lz_1[i]= Lz[i].first;
    Lz_2[i]= Lz[i].second;
  }
  fwrite( (const void*)Sr_1.data(), sizeof(double), Sr_1.size(), f );
  fwrite( (const void*)Sr_2.data(), sizeof(int), Sr_2.size(), f );
  fwrite( (const void*)Sth_1.data(), sizeof(double), Sth_1.size(), f );
  fwrite( (const void*)Sth_2.data(), sizeof(int), Sth_2.size(), f );
  fwrite( (const void*)Lz_1.data(), sizeof(double), Lz_1.size(), f );
  fwrite( (const void*)Lz_2.data(), sizeof(int), Lz_2.size(), f );

  fclose(f);
  return true;
}


/*
 * Here we create a lookup table which takes in [Sr,key] and [Sth,key] and sorts them by the first value.
 * The point is that we want to find some given Sr and Sth in order to solve the adiabatic invariants; e.g. 
 * we want to match Sr in the beginning of the growth with the Sr value after growth. To do this, we need to 
 * look up the table with Sr values in a fast manner, hence we use a sorted list for this. 
 * On the other hand, we also want to find the constants of motions that solve these equations; this is where 
 * the key comes in. For each key, we store the values of E, L, Q which can be found quickly. Then, if 
 * the Sr and Sth can be found AND they have the same key, we can be confident that this is the correct solution. 
 * The reason we do NOT use a grid methdod is that for most values of Eps, Lz, Q there are no bound orbits and hence 
 * we can not calculate  the values of Sr, Sth. Using a grid method would mean that we would not fill most of hte
 * grid space.
 * 
 */
int main(int argc, char** argv) {

  int rank;
//  MPI_Init(&argc,&argv);
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> dist(1.0, 10.0);


  // Get arguments:
  cout << "Running: ./testLookupTable <number of particles> <mass> <spin> <rMin (M)> <filename>" << endl;
  cout << "Example: ./testLookupTable 10000 500 0.998 200 lookup_0.998.dat" << endl;

  if( argc != 6 ) { 
    cout << "bad number of args" << endl; return 0;
  }
  const int nParticles = atoi(argv[1]);
  const double M = std::atof(argv[2]);
  const double spin= std::atof(argv[3]); assert( spin <= 1 );
  const double rMinLimit = std::atof(argv[4]) * M; // Input should be units of M
  const char * filename = argv[5];

  // Create metric:
  const double a = spin*M;

  // Particle sampling:
  vector<double> EpsVector(nParticles);
  vector<double> LzVector(nParticles);
  vector<double> QVector(nParticles);
  vector<pair<double,int>> Sr (nParticles);
  vector<pair<double,int>> Sth(nParticles);
  vector<pair<double,int>> LzV(nParticles);
  for( int i = 0; i < nParticles; ++i ) {
    EpsVector[i] = 0;
    LzVector[i]  = 0;
    QVector[i]   = 0;
  }
  
  // Keep a track of files:
#pragma omp parallel for schedule(dynamic,16)
  for( int i = 0; i < nParticles; ++i ) {
    if( i%1000 == 0 ) cout << "Sampling particle " << i << endl;
    double rMin, rMax, thMin, thMax;
    double consts[4]={0};
    bool success = sampleDistribution(consts, mt, M, a, rMinLimit);
    assert( success );
    const double Eps = consts[0], Lz = consts[1], Q = consts[2];
    success = helper::kerr_r_limits(Eps, Lz, Q, a, M, rMin, rMax);
    assert( success );
    success = helper::kerr_th_limits(Eps, Lz, Q, a, M, thMin, thMax);
    assert( success );
    Sr[i].first = helper::I_r(Eps, Lz, Q, a, M, rMin, rMax);
    assert( Sr[i].first == Sr[i].first ); // Check for nans
    Sr[i].second = i;
    Sth[i].first = helper::I_th(Eps, Q, Lz, a, M, rMin, rMax);
    //cout << Sth[i].first << endl;
    assert( Sth[i].first == Sth[i].first ); // Check for nans
    Sth[i].second = i;
    EpsVector[i] = Eps;
    LzVector[i]  = Lz;
    LzV[i].first  = Lz;
    LzV[i].second = i;
    QVector[i]   = Q;
  }
  // Sort the vectors
  sort(Sr.begin(), Sr.end());
  sort(Sth.begin(), Sth.end());
  sort(LzV.begin(), LzV.end());

  // Calculate the total number of particles in the whole simulation domain, including the particles we did not sample which were outside the r(MMax)
  {
    // Save results
    save_results(EpsVector,LzVector,QVector, Sr, Sth, LzV, nParticles, rMinLimit, a, M, filename);
  }

  /*print out final results into profile_0.txt*/

//  MPI_Finalize();

  return 0;
}














