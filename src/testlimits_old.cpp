// Gyoto
#include <cassert>
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
#include <random>
#include <cstdlib>
#include <cstdio>
#include <string>
//#include "mpi.h"
#include "omp.h"
#include "sampling.h"
#include <iostream>
#include <fstream>
#include "misc.h"

using namespace std;

const int starMaxIterations = 1000000;

void uniformSampling(double x[4], double u[4], double v[3], mt19937_64 & gen, const double a, const double M, const double rMinLimit, FILE * file) {
  uniform_real_distribution<> rand(0.,1.);
  // Sample E, L, Q uniformly:
  double E,L,Lz,Q,rmin,rmax;
  while( true ) {
    E = 0.7 + 0.3*rand(gen);
    L = 2.*rand(gen)*70.-70.;
    Lz= 2.*rand(gen)*70.-70.;
    //Lz= L*(rand(gen)*2.-1);
    if( abs(Lz) > abs(L) ) continue;
    Q = L*L-Lz*Lz;
    if(helper::kerr_r_limits(E, Lz, Q, a, M, rmin, rmax) == true// && rmin <= rMinLimit
                                                                 && E < 0.999 ) break; // Found bound orbit
  }
  assert( helper::kerr_r_limits(E, Lz, Q, a, M, rmin, rmax) == true );
  assert( file );
  fprintf(file,"%.16e %.16e %.16e %.16e %.16e\n",E, Lz, Q, rmin, rmax);
  helper::constants_to_pos(E,Lz,Q,a,M,x,u);
  return;
}


int main(int argc, char ** argv) {
  // Mersenne twister
  random_device rd;
  mt19937_64 gen(11);
  double x[4] = {0};
  // Test the solver:
  // Get arguments:
  cout << "Running: ./a.out <number of particles> <integrator> <M1> <spin> <rbins> <thbins> <phBins> <rMax> <filename> <adaptive 0 or 1>" << endl;
  cout << "Example: ./a.out 10000 runge_kutta_fehlberg78 1 0 1500 20 1 100 smallcore.dat 1" << endl;

  if( argc != 11 ) {
    cout << "bad number of args" << endl; return 0;
  }
  const int nParticles = atoi(argv[1]);
  const char * integrator = argv[2];
  const double M_1 = std::atof(argv[3]);
  const double spin= std::atof(argv[4]); assert( spin <= 1 );
  const int rBins = atoi(argv[5]);
  const int thBins= atoi(argv[6]);
  const int phBins= atoi(argv[7]);
  const double rMax = std::atof(argv[8]);
  const char * filename = argv[9];
  const bool adaptive = atoi(argv[10]);
  // Create metric:
  const double a_1 = spin*M_1;
  //assert( a_1 == 0 );
  // Create a grid for sampling:
  const double Rs = 2.;
  const double epsRMin = 0.00000001;
  const double rMin = (2.+sqrt(pow(Rs,2)-4*pow(spin,2)) + epsRMin)/2.;
  // Add settings for integration time
  const double integration_time = 100000;
  FILE * file = fopen(filename,"w");
  assert( file );
  for( int i = 0; i < nParticles; ++i ) {
    double x[4] = {0};
    double v[3] = {0};
    double u[4] = {0};
    // Sample particles uniformly
    uniformSampling(x,u,v,gen,a_1,M_1, rMax, file);
  }
  fclose(file);
  return 0;
}

