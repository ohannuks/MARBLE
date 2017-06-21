// Gyoto
#include "GyotoDefs.h"
#include "GyotoFactory.h"
#include "GyotoUtils.h"
#include "GyotoRegister.h"
#include "GyotoStar.h"
#include "GyotoKerrBL.h"
#include "cstdlib"
#include <iostream>
#include <chrono> 
#include <ctime>
#include <cstdio>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>
#include "helperFunctions.h"
#include <cmath>
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
#include <fitsio.h> 
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
#include "grid.h"



using namespace std;
using namespace Gyoto;

const int starMaxIterations = 10000000;


const double pi = M_PI;


// Note: I know this is very dumb, but this is only used in testconstantcore and hence the file read/write is hardcoded. see processConstantCoreData.py in the Tests folder
static
void read_results( vector<double> & EpsVector,
                   vector<double> & LzVector,
                   vector<double> & QVector,
                   vector<double> & Ndist, // This is processed by a python script; see the processData.py
                   int & Ntotal,
                   int & Nin,
                   double & a_2,
                   double & M_2,
                   const string & filename) {
  EpsVector.reserve(20000);
  LzVector.reserve(20000);
  QVector.reserve(20000);
  Ndist.reserve(20000);
  const int width = 8;
  double variables[width];
  assert( EpsVector.size() == QVector.size() );
  FILE * f = fopen(filename.c_str(), "r");
  // Scan the data
  int iter = 1;
  while( !feof(f) ) {
    for( int j = 0; j < width; ++j ) {
      double variable;
      if(fscanf(f, "%lf", &variable) == EOF) assert(1==1);
      variables[j] = variable;
    }
    if( feof(f) ) continue;
    EpsVector.push_back(variables[0]);
    LzVector.push_back(variables[1]);
    QVector.push_back(variables[2]);
    Ndist.push_back(variables[3]);
    M_2 = variables[4];
    a_2 = variables[5];
    Ntotal= variables[6];
    Nin   = variables[7];
    //cout << "Var: " << variables[0] << " " << variables[1] << " " << variables[2] << " " << variables[3] << " " << variables[4] << " " << variables[5] << " " << variables[6] << " " << variables[7] << endl;
    iter++;
  }
  cout << "M_2 a_2: " <<  M_2 << " " << a_2 << endl;
  assert( M_2 == 1 ); // Assuming these are in Gyoto units
  assert( a_2 <= 1 );
  fclose(f);
}

static
void getCoords( Gyoto::Astrobj::Star & star,
                std::vector<double> & t,
                std::vector<double> & r,
                std::vector<double> & th,
                std::vector<double> & ph,
                std::vector<double> & v_t,
                std::vector<double> & v_r,
                std::vector<double> & v_th,
                std::vector<double> & v_ph ) {
  const int nElements = star.get_nelements();
  assert( t.size() == nElements-1 );
  // Note: We use nElements - 1 because if we use the sampling procedure given by Gyoto and use constant stepping, the final value is usually t=0, which makes no sense
  for( int i = 0; i < nElements-1; ++i ) {
    double dest[8];
    star.getCoord(i, dest);
    t[i]   = dest[0];
    r[i]   = dest[1];
    th[i]  = dest[2];
    ph[i]  = dest[3];
    v_t[i] = dest[4];
    v_r[i] = dest[5];
    v_th[i]= dest[6];
    v_ph[i]= dest[7];
    if( i > 0 && t[i] < t[i-1] ) std::cerr << "ERROR: Bad t value: " << t[i] << " " << t[i-1] << " " << std::endl;
  }
  return;
}

// Integrate timelike geodesics
int integration( const Real integration_time, Gyoto::Astrobj::Star & star, SmartPointer<Gyoto::Metric::Generic> gg, 
                 Grid3D & spatialGrid, const bool saveToGrid, const Real weight ) {
  // Integrate (note: stores the coordinates; if necessary we need to optimize this)
  cout << "Integrating" << integration_time <<  endl;
  const bool newtonian = false;
  const int xFillCode = SUCCESSFUL; // TODO: Fix this, I dont actually know if I am breaking the code by doing this
  cout << "Starting integration.. " << endl;
  star.xFill( integration_time );
  cout << "Done integrating" << endl;

  if( xFillCode != SUCCESSFUL ) return UNSUCCESSFUL;

  if( saveToGrid ) {
    // Note: We use nElements-1 because for constant stepping, we have 
    const int nElements = star.get_nelements()-1;
    if( nElements < 6 ) { cout << "ERROR: Too few points for trajectory!"; return UNSUCCESSFUL; }// Insufficient points
    // Get the coordinates:
    vector<double> t(nElements); vector<double> r(nElements); vector<double> th(nElements); vector<double> ph(nElements);
    vector<double> v_t(nElements); vector<double> v_r(nElements); vector<double> v_th(nElements); vector<double> v_ph(nElements);
    getCoords(star,t,r,th,ph,v_t,v_r,v_th,v_ph);
    
    // Check that time is increasing
    for( int i = 1; i < nElements; ++i ) { if( t[i-1] > t[i] ) {cout << "ERROR bad t value: " << t[i-1] << " " << t[i] << " " << i << endl; return UNSUCCESSFUL; } }
    for( int i = 0; i < nElements; ++i ) { if( t[i] != t[i] || r[i] != r[i] || th[i] != th[i] || ph[i] != ph[i] || v_t[i] != v_t[i] || v_r[i] != v_r[i] || v_th[i] != v_th[i] || v_ph[i] != v_ph[i] ) exit(1);}
    // Make sure the particle is not traveling faster than the speed of light:
    for( int i = 1; i < nElements; ++i ) {
      Real fakePos[4] = {0};
      fakePos[0] = t[i];
      fakePos[1] = r[i];
      fakePos[2] = th[i];
      fakePos[3] = ph[i];
      Real fakeVel[4] = {0};
      fakeVel[0] = v_t[i]; 
      fakeVel[1] = v_r[i]; 
      fakeVel[2] = v_th[i];
      fakeVel[3] = v_ph[i];
      double g[4][4] = {0};
      static_cast<SmartPointer<Metric::KerrBL> >(gg) -> gmunu(g, fakePos); // Get g
      //double tetrad[4][4] = {0};
      //static_cast<SmartPointer<Metric::KerrBL> >(gg) -> ZAMOTetrad(fakePos, tetrad, true);
      //double localVel[4] = {0};
      //for( int m = 0; m < 4; ++m ) for( int mu = 0; mu < 4; ++mu ) {
      //  localVel[m] += tetrad[m][mu]*vel[mu];
      //}
      double velNorm=0;
      for( int mu = 0; mu < 4; ++mu ) for( int nu = 0; nu < 4; ++nu ) velNorm += g[mu][nu]*fakeVel[mu]*fakeVel[nu];
      if( abs(velNorm+1.) > 0.0001 ) { cout << "Bad velnorm: " << abs(velNorm+1.) << endl; return UNSUCCESSFUL;}
    }
    // Get cell properties
    const Real cell_r_width = (spatialGrid.rMax - spatialGrid.rMin)/spatialGrid.rDimensions;
    const Real cell_th_width = (spatialGrid.thMax - spatialGrid.thMin)/spatialGrid.thDimensions;
    const Real cell_ph_widph = (spatialGrid.phMax - spatialGrid.phMin)/spatialGrid.phDimensions;
    // Get number of elements
    if( nElements <= 1 ) return UNSUCCESSFUL;

    assert( t.size() > 0 );

    // Interpolate between every t units (dont care about performance for now):
    const Real tmax = t[nElements-1]; const Real tmin = t[0];
    const Real T = tmax - tmin;
    const Real R0 = spatialGrid.rMax;
    vector<double> tValues, rValues, thValues, phValues, weights, vtValues, vrValues, vthValues, vphValues, timeValues;
    vector<double> Vx, Vy, Vz, Vt;
    const int nIterations = nElements;
    tValues.reserve(nIterations);
    rValues.reserve(nIterations); thValues.reserve(nIterations); phValues.reserve(nIterations); 
    vtValues.reserve(nIterations); vrValues.reserve(nIterations); vthValues.reserve(nIterations); vphValues.reserve(nIterations);
    weights.reserve(nIterations); timeValues.reserve(nIterations);
    Vx.reserve(nIterations); Vy.reserve(nIterations); Vz.reserve(nIterations), Vt.reserve(nIterations);
    for( int i = 1; i < nElements-1; ++i ) {
      // Hack performance boost: drop out all the elements which are not in range of sampling:
      if( r[i] > spatialGrid.rMax ) { continue; }
      
      // Is in range, so sample:
      const Real t0 = t[i-1];
      const Real t1 = t[i+1];
      const Real dt = (t1-t0)/2.;
      //(t1 - t0);
      
      // NOTE: here we assume that phi does not contribute to the metric!

      // Add a weight to the spatial grid
      {
        //const double grr =   static_cast<SmartPointer<Metric::KerrBL> >(gg) -> gmunu(pos,1,1);
        //const double gthth = static_cast<SmartPointer<Metric::KerrBL> >(gg) -> gmunu(pos,2,2);
        //const double gpp =   static_cast<SmartPointer<Metric::KerrBL> >(gg) -> gmunu(pos,3,3);
        //const double volume_weight = sqrt(grr * gthth * (gpp*2.*pi));
        tValues.push_back(  t[i] );
        rValues.push_back(  r[i] );
        thValues.push_back( th[i]);
        phValues.push_back( ph[i] );
        vtValues.push_back(  v_t[i] );
        vrValues.push_back(  v_r[i] );
        vthValues.push_back( v_th[i] );
        vphValues.push_back( v_ph[i] );
        // Add velocity values
        timeValues.push_back(dt);
        weights.push_back(dt * weight); // Additional weight based on particle distribution
      }
    }

    // Input values
    assert( rValues.size() == thValues.size() && rValues.size() == phValues.size() );
    #pragma omp critical
    {
      const double M = 1;
      const double a = static_cast<SmartPointer<Metric::KerrBL> >(gg) -> spin();
      const int sampleSuccessful = spatialGrid.sample(tValues, rValues, thValues, phValues, vtValues, vrValues, vthValues, vphValues, weights, a, M, T);
//      for( int i = 0; i < rValues.size(); ++i ) {
//        if( rValues[i]  >= 0 ) {
//          const int sampleSuccessful = spatialGrid.sample(0, rValues[i], thValues[i], phValues[i], vtValues[i], vrValues[i], vthValues[i], vphValues[i], timeValues[i], a, M);
//        }
//      }
    }
  }
  return SUCCESSFUL;
}


int main(int argc, char** argv) {

  int rank;
//  MPI_Init(&argc,&argv);
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Get arguments:
  cout << "Running: ./rungyoto <integrator> <rbins> <thbins> <phBins> <rMax> <filename> <adaptive 0 or 1> <filename in>" << endl;
  cout << "Example: ./rungyoto runge_kutta_fehlberg78 200 30 1 200 Data/degeneracygrid 1 Data/degeneracydistributiongyoto.dat" << endl;

  if( argc < 8 ) { 
    cout << "bad number of args" << endl; return 0;
  }
  const char * integrator = argv[1];
  const int rBins = atoi(argv[2]);
  const int thBins= atoi(argv[3]);
  const int phBins= atoi(argv[4]);
  const double rMax = std::atof(argv[5]);
  const char * filename = argv[6];
  const bool adaptive = atoi(argv[7]);
  const char * fin      = argv[8];
  
  vector<double> EpsVector, LzVector, QVector, Ndist;
  int Ntotal=0, Nin=0;
  double a_2, M_2;
  read_results(EpsVector, LzVector, QVector, Ndist, Ntotal, Nin, a_2, M_2, fin);
  const int nParticles = EpsVector.size();

  assert( a_2 < 1 && M_2 == 1);
  const double spin = a_2;
  // Create a grid for sampling:
  const double Rs = 2.;
  const double epsRMin = 0;
  const double rMin = (2.+sqrt(pow(Rs,2)-4*pow(spin,2)) + epsRMin)/2.;
  const double dr   = (rMax - rMin)/(double)rBins;
  Grid3D spatialGrid( rBins, thBins, phBins, rMax, rMin );

  int skippedParticles = 0;
  const double integration_time = 100000;

  // Main sampling
#pragma omp parallel for schedule(dynamic,16) reduction(+:skippedParticles)
  for( int i = 0; i < nParticles; ++i ) {
    if( i%1000 == 0 ) cout << "SAMPLING PARTICLE " << i << endl;
    double x[4] = {0};
    double v[3] = {0};
    double u[4] = {0};
    double consts[4]={0};
    double rmin, rmax, thmin, thmax;
    // Make sure this is a bound orbit:
    const double Eps = EpsVector[i];  
    const double L   = LzVector[i]; 
    const double Q   = QVector[i]; 
    bool success = (helper::kerr_r_limits(Eps, L, Q, a_2, M_2, rmin,rmax) && helper::kerr_th_limits(Eps, L, Q, a_2, M_2, thmin,thmax));
    assert( success );
    // Sample the particle:
    // Put the particle in the correct position:
    int dir = 2*(rand()%2)-1;
    success = helper::constants_to_pos(Eps,L, Q, a_2, M_2, x, u, dir);
    assert( success );
    v[0] = u[1]/u[0];
    v[1] = u[2]/u[0];
    v[2] = u[3]/u[0];
    if( i%1000 == 0 ) cout << "Setting particle with Eps, Lz, Q, x, u: " << Eps << " " << L << " " << Q << " " << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << " " << u[0] << " " << u[1] << " " << u[2] << " " << u[3] << endl;
    // Sample the particle!
    // Create a particle:
    SmartPointer<Gyoto::Metric::Generic> gg (new Gyoto::Metric::KerrBL());
    static_cast<SmartPointer<Metric::KerrBL> >(gg) -> spin(spin);
    static_cast<SmartPointer<Metric::KerrBL> >(gg) -> mass(1.);
    Gyoto::Astrobj::Star darkMatterParticle( static_cast<SmartPointer<Metric::Generic> >(gg),
         0.001,
         x, v);
    {
      darkMatterParticle.integrator(integrator);
      if( adaptive ) {darkMatterParticle.adaptive(true);  darkMatterParticle.deltaMaxOverR(0.08); darkMatterParticle.deltaMax(dr*0.5); } // set adaptive integration
      else           {darkMatterParticle.adaptive(false); darkMatterParticle._delta(0.05); }// No adaptive stepping
      darkMatterParticle.maxiter(starMaxIterations); // Max number of iterations -> a lot
    }
    const double weight = Ndist[i];
    int integrationSuccessful = integration(integration_time, darkMatterParticle, gg, spatialGrid, true, weight);
    if( integrationSuccessful == SUCCESSFUL ) { if(i%1000 == 0) cout << "INTEGRATION SUCCESSFUL" << endl; }
    if( integrationSuccessful == UNSUCCESSFUL ) { 
      cout << "INTEGRATION UNSUCCESSFUL" << endl;
      ++skippedParticles;
    }
  }

  cout << "Done integrating, saving results.." << endl;
  // Calculate number of unsampled particles:
  {
    // Normalize grid 
    //spatialGrid.normalize_grid(Ntotal);
    // Save
    string filename2 = filename;
    cout << "Saving to grid.." << endl;
    spatialGrid.save_grid(filename2);
  }
  cout << "Skipped particles: " << skippedParticles << " Out of " << nParticles << endl;

  /*print out final results into profile_0.txt*/

//  MPI_Finalize();

  return 0;
}














