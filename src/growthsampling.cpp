// Gyoto
#include "cstdlib"
#include <iostream>
#include <chrono> 
#include <ctime>
#include <cstdio>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>
#include "helperFunctions.h"
#include <cmath>
#include <cassert>
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
#include "misc.h"
#include "haloFunctions.h"

using namespace std;
const int starMaxIterations = 1000000;

const double pi = M_PI;




int main(int argc, char** argv) {

  int rank;
//  MPI_Init(&argc,&argv);
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Get arguments:
  cout << "Running: ./growthsampling <number of particles> <M1> <M2> <spin> <rMax> <filename> <maxwellian #x> <#y> <#z> <maxwellian dispersion> <density profile>" << endl;
  cout << "Example: ./growthsampling 10000 1 2 1 150 griddata.dat 0 0 0 0.0006671 2" << endl;

  if( argc < 12 ) { 
    cout << "bad number of args" << endl; return 0;
  }
  const int nParticles = atoi(argv[1]);
  const double M_1 = std::atof(argv[2]);
  const double M_2 = std::atof(argv[3]);
  const double spin= std::atof(argv[4]); assert( spin <= 1 );
  const double rMax = std::atof(argv[5]);
  const char * filename = argv[6];
  const double mean_v[3] = {std::atof(argv[7]),std::atof(argv[8]),std::atof(argv[9])}; // Mean velocity of the Maxwellian distribution
  // TODO: This is actually no longer true, r_influence should be a function of M_1
  const double sigma = std::atof(argv[10]);
  const int densityProfileN = atoi(argv[11]);

  sampling::Profile profile;
  switch( densityProfileN ) {
    case 0:
      profile = sampling::Profile::constant;
      break;
    case 1:
      profile = sampling::Profile::r1;
      break;
    case 2:
      profile = sampling::Profile::r2;
      break;
    case 3:
      profile = sampling::Profile::r3;
      break;
    default:
      assert(1==0);
      break;
  }
  const double rMinInfluence = M_1 / pow(sigma,2);
  const double rMaxInfluence = M_2 / pow(sigma,2);

  // Mersenne twister
  random_device rd;
  mt19937_64 gen(11);
  // Normal distribution c++11
  const double mean = 0;
  normal_distribution<double> normal_dist(mean, sigma);

  // Create metric:
  //TODO: This is probably incorrect
  const double a_1 = spin*M_1;
  const double a_2 = spin*M_2;

  // Create a grid for sampling:
  const double Rs = 2.;
  const double epsRMin = 0;
  const double rMin = (2.+sqrt(pow(Rs,2)-4*pow(spin,2)) + epsRMin)/2.;

  // Keep track of the sampling process:
  // All of the particles that get sampled within r<rMaxSample
  long long Nsampled = 0; 
  // All of the particles that get sampled AND survive adiabatic growth
  long long Nin      = 0; 
  
  // Particle sampling:
  vector<double> EpsVectorNew(nParticles);
  vector<double> LzVectorNew(nParticles);
  vector<double> QVectorNew(nParticles);
  vector<double> EpsVector(nParticles);
  vector<double> LzVector(nParticles);
  vector<double> QVector(nParticles);
  vector<double> MStartVector(nParticles);
  vector<double> aStartVector(nParticles);
  for( int i = 0; i < nParticles; ++i ) {
    EpsVector[i] = 0;
    LzVector[i]  = 0;
    QVector[i]   = 0;
    EpsVectorNew[i] = 0;
    LzVectorNew[i]  = 0;
    QVectorNew[i]   = 0;
    MStartVector[i] = 0;
    aStartVector[i] = 0;
  }
  
  // Keep a track of files:
  const double MMax = 5;
  const double rMaxSample = MMax/pow(sigma,2);
#pragma omp parallel for schedule(dynamic,16) reduction(+ : Nsampled, Nin )
  for( int i = 0; i < nParticles; ++i ) {
    cout << "SAMPLING PARTICLE " << i << endl;
    double x[4] = {0};
    double v[3] = {0};
    double u[4] = {0};
    double consts[4]={0};
    double consts_1[4] = {0};
    double consts_2[4] = {0};

    // Do adiabatic growth
    double rmin, rmax, thmin,thmax, rMinOld, rMaxOld, M_start=0, a_start=0;
    while( true ) {
      while( true ) {
        // R sampling
        sampling::r_sampling(rMinInfluence, rMaxInfluence,gen, x, profile,rMaxSample);
        assert( x[1]*pow(sigma,2) < MMax );
        M_start = x[1]*pow(sigma,2);
        a_start = spin*M_start;
        // Velocity sampling
        sampling::maxwellian_sampling(normal_dist, gen, mean_v, x, v);
        helper::findVt(x, v, a_1, M_1, u); // Normalize u^mu

        // Check to make sure the orbit is bound:
        {
          helper::constants(x,u,a_start,M_start,consts);
          const double Eps = consts[0];
          const double L   = consts[1];
          const double Q   = consts[2];
          if( Eps > 1 || helper::kerr_r_limits(Eps, L, Q,a_start, M_start, rmin,rmax) == false ) {
            Nsampled=Nsampled+1;
            continue;
          }
        }
        break;
      }
      for( int i = 0; i < 4; ++i ) consts_1[i] = consts[i];
      rMinOld = rmin; rMaxOld = rmax;

      // Grow black hole
      double result = 0;
      Solver solver = Solver::hybrid;
      int growthSuccessful = helper::grow_black_hole(consts, a_start, a_2, M_start, M_2, solver, consts_2, result);
      // Check if the BH growth had a solution:
      if( growthSuccessful == SUCCESSFUL ) {
        // Get orbit:
        {
          const double Eps = consts_2[0];
          const double L = consts_2[1];
          const double Q = consts_2[2];
          if( Eps > 1 || helper::kerr_r_limits(Eps, L, Q, a_2, M_2, rmin,rmax) == false ) {
            assert(1==1); // This should never happen
            continue;
          }
        }
        // Make sure the particle is sampled on the grid
        if( rmin/M_2 > rMax ) {
          continue;
        }
        // Particle sampling successful!
        break;
      }
    }
    // Make sure this is a bound orbit:
    const double Eps = consts_2[0];
    const double L = consts_2[1];
    const double Q = consts_2[2];
    assert(helper::kerr_r_limits(Eps, L, Q, a_2, M_2, rmin,rmax) == true && helper::kerr_th_limits(Eps, L, Q, a_2, M_2, thmin,thmax));
    assert(helper::kerr_r_limits(consts_1[0], consts_1[1], consts_1[2], a_start, M_start, rmin,rmax) == true);
    // Sample the particle
    EpsVectorNew[i] = consts_2[0];
    LzVectorNew[i]  = consts_2[1];
    QVectorNew[i]   = consts_2[2];
    EpsVector[i]    = consts_1[0];
    LzVector[i]     = consts_1[1];
    QVector[i]      = consts_1[2];
    MStartVector[i] = M_start;
    aStartVector[i] = a_start;
    Nin             = Nin+1;
    cout << "Inserted rLimits and rLimitsOld: " << rmin << " " << rmax << " " << thmin << " " << thmax << " " << rMinOld << " " << rMaxOld << " " << consts[0] << " " << consts_2[0] << " " << M_start <<  " " << endl;
  }
  // Calculate the total number of particles in the whole simulation domain, including the particles we did not sample which were outside the r(MMax)
  {
    const int n = densityProfileN;
    const long double r1    = M_1/pow(sigma,2);
    const long double rX    = MMax/pow(sigma,2);
    const long double r2    = M_2/pow(sigma,2);
    const long double Ntotal= (long double)Nsampled * dM(n,r1,r2) / dM(n, r1, rX);
    cout << "NTotal: " <<  Ntotal << " Nsampled: " << Nsampled <<  " dM(r1,r2) " << dM(n,r1,r2) << " dM(r1,rX) " << dM(n,r1,rX) << endl;
    const string filename2 = filename;
    save_results(EpsVector, LzVector, QVector, 
                 EpsVectorNew, LzVectorNew, QVectorNew, 
                 MStartVector, aStartVector,
                 Ntotal, Nin, a_1, M_1, a_2, M_2, filename2); // For normalization; density goes as numberOfParticles // FIX! See the TODO above!!!!!! This is actually important
  }

  /*print out final results into profile_0.txt*/

//  MPI_Finalize();

  return 0;
}














