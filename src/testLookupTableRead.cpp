#include <cassert>
#include "cstdlib"
#include <algorithm>
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
#include "haloFunctions.h"
#include "misc.h"
#include "lookupTable.h"
#include "GrowBlackHole.h"

using namespace std;
const int starMaxIterations = 1000000;

const double pi = M_PI;



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
  


  // Get arguments:
  cout << "Running: ./growthsampling <number of particles> <M1> <M2> <spin> <rMax> <filename> <maxwellian #x> <#y> <#z> <maxwellian dispersion> <density profile> <lookuptable>" << endl;
  cout << "Example: ./growthsampling 10000 1 2 0.998 150 griddata.dat 0 0 0 0.0006671 2 lookup_0.998.dat" << endl;

  if( argc != 13 ) { 
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
  const char * lookupfname  = argv[12];

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
  // Create metric:
  const double a_1 = spin*M_1;
  const double a_2 = spin*M_2;

  // Mersenne twister
  random_device rd;
  mt19937_64 gen(11);
  // Normal distribution c++11
  const double mean = 0;
  normal_distribution<double> normal_dist(mean, sigma);
  
  // Read the lookup table
  LookupTable lookup;
  read_lookup(lookup, lookupfname );
  
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

  // Keep track of the sampling process:
  // All of the particles that get sampled within r<rMaxSample
  int Nsampled = 0; 
  // All of the particles that get sampled AND survive adiabatic growth
  int Nin      = 0; 
  const double MMax = 4;
  const double rMaxSample = MMax/pow(sigma,2);

  GrowBlackHole( nParticles,
                 spin,
                 a_1,
                 M_1,
                 a_2,
                 M_2,
                 rMinInfluence,
                 rMaxInfluence,
                 rMaxSample,
                 MMax,
                 sigma,
                 mean_v,
                 rMax,
                 profile,
                 gen,
                 normal_dist,
                 lookup,
                 Nsampled,
                 Nin,
                 EpsVectorNew,
                 LzVectorNew,
                 QVectorNew,
                 EpsVector,
                 LzVector,
                 QVector,
                 MStartVector,
                 aStartVector);

  // Calculate the total number of particles in the whole simulation domain, including the particles we did not sample which were outside the r(MMax)
  {
    const int n = densityProfileN;
    const double r1    = M_1/pow(sigma,2);
    const double rX    = MMax/pow(sigma,2);
    const double r2    = M_2/pow(sigma,2);
    const double Ntotal= Nsampled * dM(n,r1,r2) / dM(n, r1, rX);
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














