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
#include <cstdlib>
#include <cstdio>
#include <string>
//#include "mpi.h"
#include "omp.h"
#include "sampling.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>
#include <vector>
#include <utility>

using namespace std;

// Global variables
vector<double> EpsVector;
vector<double> LzVector;
vector<double> QVector;
vector<pair<double,int>> Sr ;
vector<pair<double,int>> Sth;
vector<pair<double,int>> LzV;
int i;


const int starMaxIterations = 1000000;

const double pi = M_PI;

double
g (double *k, size_t dim, void *params)
{
  assert( dim == 3 );
  double * paramsDouble = (double*)params;
  const double Eps = k[0], Lz = k[1], Q = k[2], K=k[3];
  const double M         = paramsDouble[0];
  const double a         = paramsDouble[1];
  const double rMinLimit = paramsDouble[2];
  double rMin, rMax, thMin, thMax;
  if( helper::kerr_r_limits(Eps, Lz, Q, a, M, rMin, rMax) && helper::kerr_th_limits(Eps, Lz, Q, a, M, thMin, thMax) && rMin < rMinLimit ) { 
    double _Sr = helper::I_r(Eps, Lz, Q, a, M, rMin, rMax),
           _Sth= helper::I_th(Eps, Q, Lz, a, M, thMin, thMax),
           _Lz = Lz;
    Sr.push_back(  make_pair(_Sr,  i));
    Sth.push_back( make_pair(_Sth, i));
    LzV.push_back( make_pair(_Lz,  i));
    EpsVector.push_back(Eps);
    LzVector.push_back(Lz);
    QVector.push_back(Q);
    i += 1;
    //0.990848 0.00561543 0.00353704
    double total = _Sr/0.995069 + _Sth/0.00561543 + _Lz/0.00353704;
    return total;
  } 
  return 0;
}



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
                  const double rMinLimit,
                  const double a,
                  const double M,
                  const char * filename
                 ) {
  const int nParticles = EpsVector.size();
  assert( EpsVector.size() == nParticles );
  assert( EpsVector.size() == Sth.size() );
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


void
display_results (char *title, double result, double error)
{
  printf ("%s ==================\n", title);
  printf ("result = % .6f\n", result);
  printf ("sigma  = % .6f\n", error);
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
  EpsVector.reserve(nParticles);
  LzVector.reserve(nParticles);
  QVector.reserve(nParticles);
  Sr .reserve(nParticles);
  Sth.reserve(nParticles);
  LzV.reserve(nParticles);
  for( int i = 0; i < nParticles; ++i ) {
    EpsVector[i] = 0;
    LzVector[i]  = 0;
    QVector[i]   = 0;
  }
  
  // MONTE CARLO IMPORTANCE SAMPLING
  {
    double res, err;
    const double Emin  = 0.9,    Emax  = 1., 
                 Lzmin = -50.*M, Lzmax = 50.*M, 
                 Qmin  = 0,      Qmax  = Lzmax*Lzmax;
    double xl[3] = { Emin, Lzmin, Qmin };
    double xu[3] = { Emax, Lzmax, Qmax };
    // Initializing everything for gsl
    const gsl_rng_type *T;
    gsl_rng *r;
    const int ndim = 3; // Eps, Lz, Q
    double params[3] = {M, a, rMinLimit}; // See function g
    gsl_monte_function G = { &g, ndim, &params[0] };
    size_t calls = nParticles;
    gsl_rng_env_setup ();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    // Monte carlo sampling
    {
      gsl_monte_miser_state *s = gsl_monte_miser_alloc (ndim);
      gsl_monte_miser_integrate (&G, xl, xu, ndim, calls, r, s,
                                 &res, &err);
      gsl_monte_miser_free (s);
      display_results ("miser", res, err);
    }
    gsl_rng_free (r);
  }

  // Sort the vectors
  sort(Sr.begin(), Sr.end());
  sort(Sth.begin(), Sth.end());
  sort(LzV.begin(), LzV.end());
  const double SrMax = Sr[Sr.size()-1].first, SthMax = Sth[Sth.size()-1].first, LzMax = LzV[LzV.size()-1].first,
               SrMin = Sr[0].first,           SthMin = Sth[0].first,            LzMin = LzV[0].first;
  const double total = (SrMax - SrMin) + (SthMax - SthMin) + (LzMax - LzMin);

  // Calculate the total number of particles in the whole simulation domain, including the particles we did not sample which were outside the r(MMax)
  {
    // Save results
    save_results(EpsVector,LzVector,QVector, Sr, Sth, LzV, rMinLimit, a, M, filename);
  }
  cout << "Total " << (SrMax-SrMin)/total << " " << (SthMax - SthMin) / total << " " << (LzMax-LzMin) / total << endl;
  cout << "nParticles: " << EpsVector.size() << endl;

  /*print out final results into profile_0.txt*/

//  MPI_Finalize();

  return 0;
}














