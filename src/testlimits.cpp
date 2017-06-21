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
// FITS I/O
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
#include "sampling.h"
#include <cstdlib>
#include <cstdio>
#include <string>
//#include "mpi.h"
#include "omp.h"
#include <iostream>
#include <fstream>
#include "misc.h"

using namespace std;

const int starMaxIterations = 1000000;

void uniformSampling(double x[4], double u[4], double v[3], mt19937_64 & gen, const double a, const double M, const double rMinLimit, 
                     vector <double > &EpsVector, 
                     vector <double > &LzVector, 
                     vector <double > &QVector, 
                     vector <double > &rMinVector, 
                     vector <double > &rMaxVector, 
                     vector <double > &thMinVector,
                     vector <double > &thMaxVector, 
                     const int i) {
  uniform_real_distribution<> rand(0.,1.);
  // Sample E, L, Q uniformly:
  double E,L,Lz,Q,rmin,rmax, thmin,thmax;
  while( true ) {
    E = 0.7 + 0.3*rand(gen);
    L = 2.*rand(gen)*70.-70.;
    Lz= 2.*rand(gen)*70.-70.;
    //Lz= L*(rand(gen)*2.-1);
    if( abs(Lz) > abs(L) ) continue;
    Q = L*L-Lz*Lz;
    if(helper::kerr_r_limits(E, Lz, Q, a, M, rmin, rmax) == true && helper::kerr_th_limits(E, Lz, Q, a, M, thmin, thmax) == true // && rmin <= rMinLimit
                                                                 && E < 0.999 ) break; // Found bound orbit
  }
  EpsVector[i] = E;
  LzVector[i]  = Lz;
  QVector[i]   = Q;
  rMinVector[i]= rmin;
  rMaxVector[i]= rmax;
  thMinVector[i]= thmin;
  thMaxVector[i]= thmax;
  return;
}

int main(int argc, char ** argv) {
  // Mersenne twister
  random_device rd;
  mt19937_64 gen(11);
  double x[4] = {0};
  // Test the solver:
  // Get arguments:
  cout << "Running: ./testlimits <number of particles> <M1> <spin> <rMax> <filename>" << endl;
  cout << "Example: ./testlimits 10000 1 0 100 limits.dat" << endl;
  cout << argv[0] << endl;
  if( argc < 6 ) {
    cout << "bad number of args" << endl; return 0;
  }
  const int nParticles = atoi(argv[1]);
  const double M_1 = std::atof(argv[2]);
  const double spin= std::atof(argv[3]); assert( spin <= 1 );
  const double rMax = std::atof(argv[4]);
  const char * filename = argv[5];

  // Create metric:
  const double a_1 = spin*M_1;

  // Create a grid for sampling:
  const double Rs = 2.;
  const double epsRMin = 0.00000001;
  const double rMin = (2.+sqrt(pow(Rs,2)-4*pow(spin,2)) + epsRMin)/2.;

  // Add settings for integration time
  vector <double > EpsVector(nParticles);
  vector <double > LzVector(nParticles);
  vector <double > QVector(nParticles); 
  vector <double > rMinVector(nParticles);
  vector <double > thMinVector(nParticles);
  vector <double > rMaxVector(nParticles); 
  vector <double > thMaxVector(nParticles);
  vector <double > MStartVector(nParticles);
  vector <double > aStartVector(nParticles);
#pragma omp parallel for schedule(dynamic,10)
  for( int i = 0; i < nParticles; ++i ) {
    double x[4] = {0};
    double v[3] = {0};
    double u[4] = {0};
    // Sample particles uniformly
    uniformSampling(x,u,v,gen,a_1,M_1, rMax, EpsVector,LzVector,QVector,rMinVector,rMaxVector,thMinVector,thMaxVector, i);
  }
  for( int i = 0; i < nParticles; ++i ) { MStartVector[i] = M_1; aStartVector[i] = a_1; }
  save_results(EpsVector, LzVector, QVector, 
               EpsVector, LzVector, QVector, 
               MStartVector, aStartVector,
               0, 0, a_1, M_1, a_1, M_1, filename); // For normalization; density goes as numberOfParticles // FIX! See the TODO above!!!!!! This is actually important
  return 0;
}

