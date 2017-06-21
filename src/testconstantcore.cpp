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
#include "grid.h"
#include "sampling.h"
#include <cstdlib>
#include <cstdio>
#include <string>
//#include "mpi.h"
#include "omp.h"
#include <sampling.h>
#include <iostream>
#include <fstream>

using namespace std;
using namespace Gyoto;

const int starMaxIterations = 100000000;

// Note: I know this is very dumb, but this is only used in testconstantcore and hence the file read/write is hardcoded. see processConstantCoreData.py in the Tests folder
static
void read_results( vector<double> & EpsVector,
                   vector<double> & LzVector,
                   vector<double> & QVector,
                   vector<double> & rMinVector,
                   vector<double> & rMaxVector,
                   vector<double> & thMinVector,
                   vector<double> & thMaxVector,
                   vector<double> & degeneraciesVector, // This is processed by a python script; see the processData.py
                   int & numberOfParticles,
                   double & a_1,
                   double & M_1,
                   const string & filename) {
  EpsVector.reserve(20000);
  LzVector.reserve(20000);
  QVector.reserve(20000);
  rMinVector.reserve(20000);
  rMaxVector.reserve(20000);
  thMinVector.reserve(20000);
  thMaxVector.reserve(20000);
  degeneraciesVector.reserve(20000);
  const int width = 8;
  assert( EpsVector.size() == QVector.size() && QVector.size() == rMinVector.size() );
  FILE * f = fopen(filename.c_str(), "r");
  double variables[width] = {0};
  for( int j = 0; j < width; ++j ) {
    double variable;
    if(fscanf(f, "%lf", &variable) == EOF) assert(1==1);
    variables[j] = variable;
  }
  // Retrieve attributes
  numberOfParticles = variables[0];
  a_1               = variables[1];
  M_1               = variables[2];
  // Scan the data
  int iter = 1;
  while( !feof(f) ) {
    for( int j = 0; j < width; ++j ) {
      double variable;
      if(fscanf(f, "%lf", &variable) == EOF) assert(1==1);
      variables[j] = variable;
    }
    EpsVector.push_back(variables[0]);
    LzVector.push_back(variables[1]);
    QVector.push_back(variables[2]);
    rMinVector.push_back(variables[3]);
    rMaxVector.push_back(variables[4]);
    thMinVector.push_back(variables[5]);
    thMaxVector.push_back(variables[6]);
    degeneraciesVector.push_back(variables[7]);
    cout << "At iteration " << iter << " " << EpsVector[EpsVector.size()-1] << " " << filename << endl;
    iter++;
  }
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
                 Grid3D & spatialGrid, const bool saveToGrid ) {
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
      }
    }

    // Input values
    assert( rValues.size() == thValues.size() && rValues.size() == phValues.size() );
    #pragma omp critical
    {
      const double M = 1;
      const double a = static_cast<SmartPointer<Metric::KerrBL> >(gg) -> spin();
      const int sampleSuccessful = spatialGrid.sample(tValues, rValues, thValues, phValues, vtValues, vrValues, vthValues, vphValues, timeValues, a, M, T);
//      for( int i = 0; i < rValues.size(); ++i ) {
//        if( rValues[i]  >= 0 ) {
//          const int sampleSuccessful = spatialGrid.sample(0, rValues[i], thValues[i], phValues[i], vtValues[i], vrValues[i], vthValues[i], vphValues[i], timeValues[i], a, M);
//        }
//      }
    }
  }
  return SUCCESSFUL;
}


void uniformSampling(double x[4], double u[4], double v[3], mt19937_64 & gen, const double a, const double M, const double rMinLimit) {
  uniform_real_distribution<> rand(0.,1.);
  // Sample E, L, Q uniformly:
  double E,L,Lz,Q,rmin,rmax;
  while( true ) {
    E = 0.9 + 0.1*rand(gen);
    L = 2.*rand(gen)*70.-70.;
    //Lz= 2.*rand(gen)*70.-70.;
    Lz= L*(rand(gen)*2.-1);
    if( abs(Lz) > abs(L) ) continue;
    Q = L*L-Lz*Lz;
    if(helper::kerr_r_limits(E, Lz, Q, a, M, rmin, rmax) == true// && rmin <= rMinLimit
                                                                 && E < 0.999 ) break; // Found bound orbit
  }
  helper::constants_to_pos(E,Lz,Q,a,M,x,u);
  v[0] = u[1]/u[0];
  v[1] = u[2]/u[0];
  v[2] = u[3]/u[0];
  return;
}


int main(int argc, char ** argv) {
  // Mersenne twister
  random_device rd;
  mt19937_64 gen(11);
  double x[4] = {0};
  // Test the solver:
  // Get arguments:
  cout << "Running: ./a.out <integrator>  <rbins> <thbins> <phBins> <rMax> <filename> <adaptive 0 or 1> <filename to read>" << endl;
  cout << "Example: ./a.out runge_kutta_fehlberg78 1500 20 1 100 smallcore.dat 1 fin.dat" << endl;

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

  // Create metric:
  double a_1, M_1;

  // Particle sampling:
  vector<double> EpsVector;
  vector<double> LzVector;
  vector<double> QVector;
  vector<double> rMinVector;
  vector<double> rMaxVector;
  vector<double> thMinVector;
  vector<double> thMaxVector;
  vector<double> degeneraciesVector;
  
  // Read data
  int numberOfSimulationParticles = 0;
  read_results(EpsVector, LzVector, QVector, rMinVector, rMaxVector, thMinVector, thMaxVector,degeneraciesVector,numberOfSimulationParticles,a_1, M_1, fin);
  const int nParticles = EpsVector.size();
  //assert( a_1 == 0. );
  assert( M_1 == 1. ); // For now this is how this is built; there is no conversion to the units M==1 as assumed by Gyoto
  
  // Create a grid for sampling:
  const Real Rs = 2.;
  const Real epsRMin = 0.00000001;
  const double spin = a_1/M_1;
  const Real rMin = (2.+sqrt(pow(Rs,2)-4*pow(spin,2)) + epsRMin)/2.;
  Grid3D spatialGrid( rBins, thBins, phBins, rMax, rMin );
  
  // Add settings for integration time
  int skippedParticles = 0;
  const double integration_time = 100000;
#pragma omp parallel for schedule(dynamic,1) reduction(+:skippedParticles)
  for( int i = 0; i < nParticles; ++i ) {
    SmartPointer<Gyoto::Metric::Generic> gg (new Gyoto::Metric::KerrBL());
    static_cast<SmartPointer<Metric::KerrBL> >(gg) -> spin(spin);
    static_cast<SmartPointer<Metric::KerrBL> >(gg) -> mass(1.);
    // Get the relevant constants 
    const double Eps = EpsVector[i], Lz = LzVector[i], Q = QVector[i];
    double x[4] = {0};
    double v[3] = {0};
    double u[4] = {0};
    // Put the particle in the correct position:
    int dir = 2*(rand()%2)-1;
    double rmin,rmax, thmin, thmax;
    if( helper::kerr_r_limits(Eps,Lz,Q,a_1,M_1,rmin,rmax) == false) {  cout << "Skipped particle " << i << endl; ++skippedParticles; continue;}
    if( helper::constants_to_pos(Eps,Lz, Q, a_1, M_1, x, u, dir) == false ) {  cout << "Skipped particle " << i << endl; ++skippedParticles; continue;}
    assert(helper::kerr_r_limits(Eps,Lz,Q,a_1,M_1,rmin,rmax));
    assert(helper::kerr_th_limits(Eps,Lz,Q,a_1,M_1,thmin,thmax));
    if(helper::constants_to_pos(Eps,Lz, Q, a_1, M_1, x, u, dir) != true ) continue;
    v[0] = u[1]/u[0];
    v[1] = u[2]/u[0];
    v[2] = u[3]/u[0];
    cout << "Setting particle with Eps, Lz, Q, x, u: " << Eps << " " << Lz << " " << Q << " " << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << " " << u[0] << " " << u[1] << " " << u[2] << " " << u[3] << endl;
    // Sample the particle!
    // Create a particle:
    Gyoto::Astrobj::Star darkMatterParticle( static_cast<SmartPointer<Metric::Generic> >(gg),
         0.001,
         x, v);
    {
      darkMatterParticle.integrator(integrator);
      if( adaptive ) {darkMatterParticle.adaptive(true);  darkMatterParticle.deltaMaxOverR(0.05); } // set adaptive integration
      else           {darkMatterParticle.adaptive(false); darkMatterParticle._delta(0.05); }// No adaptive stepping
      darkMatterParticle.maxiter(starMaxIterations); // Max number of iterations -> a lot
    }
    const double degeneracy = degeneraciesVector[i];
    int integrationSuccessful = integration(integration_time, darkMatterParticle, gg, spatialGrid, true);
    if( integrationSuccessful == SUCCESSFUL ) cout << "INTEGRATION SUCCESSFUL" << i << " " << nParticles << endl;
    if( integrationSuccessful == UNSUCCESSFUL ) { 
      cout << "INTEGRATION UNSUCCESSFUL " << i << " " << nParticles << " " << endl;
    }
  }
  cout << "Done integrating, saving results.." << endl;
  // Calculate number of unsampled particles:
  {
    // Normalize grid 
#warning no normalization
    spatialGrid.normalize_grid(1.);
    // Save
    string filename2 = filename;
    filename2.append(".dat");
    cout << "Saving to grid.." << endl;
    spatialGrid.save_grid(filename2);
  }
  cout << "Skipped particles: " << skippedParticles << " Out of " << nParticles << endl;
  cout << "Spin is: " << spin << endl;
  return 0;
}

