#include <cassert>
#include "helperFunctions.h"
#include <iostream>
#include "sampling.h"
#include "lookupTable.h"
using namespace std;

int GrowBlackHole( const int nParticles,
                   const double spin,
                   const double a_1,
                   const double M_1,
                   const double a_2,
                   const double M_2,
                   const double rMinInfluence,
                   const double rMaxInfluence,
                   const double rMaxSample,
                   const double MMax,
                   const double sigma,
                   const double mean_v[3],
                   const double rMaxLimit,
                   sampling::Profile & profile,
                   mt19937_64 & gen,
                   normal_distribution <double > & normal_dist,
                   LookupTable & lookup,
                   int & Nsampled,
                   int & Nin,
                   vector<double> & EpsVectorNew,
                   vector<double> & LzVectorNew,
                   vector<double> & QVectorNew,
                   vector<double> & EpsVector,
                   vector<double> & LzVector,
                   vector<double> & QVector,
                   vector<double> & MStartVector,
                   vector<double> & aStartVector,
                   const bool sampleUnsuccessfulGrowths
) {

  int lookupGrowthSuccesses = 0;
  int traditionalGrowthSuccesses = 0;
  
  Nsampled = 0;
  Nin = 0;
  if( EpsVectorNew.size() != nParticles ) {
    EpsVectorNew.resize(nParticles);
    LzVectorNew.resize(nParticles);
    QVectorNew.resize(nParticles);
    EpsVector.resize(nParticles);
    LzVector.resize(nParticles);
    QVector.resize(nParticles);
    MStartVector.resize(nParticles);
    aStartVector.resize(nParticles);
  }
#pragma omp parallel for schedule(dynamic,16) reduction(+ : Nsampled, Nin, lookupGrowthSuccesses, traditionalGrowthSuccesses )
  for( int i = 0; i < nParticles; ++i ) {
    cout << "SAMPLING PARTICLE " << i << endl;
    double x[4] = {0};
    double v[3] = {0};
    double u[4] = {0};
    double consts[4]={0};
    double consts_1[4] = {0};
    double consts_2[4] = {0};
    int growthSuccessful;

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
          if( Eps > 1 || helper::kerr_r_limits(Eps, L, Q,a_start, M_start, rmin,rmax) == false || helper::kerr_th_limits(Eps,L,Q,a_start,M_start,thmin,thmax) == false ) {
            Nsampled=Nsampled+1;
            continue;
          }
        } 
        break;
      }
      for( int i = 0; i < 4; ++i ) consts_1[i] = consts[i];
      rMinOld = rmin; rMaxOld = rmax;

      // Grow black hole
      const double Eps_1 = consts_1[0], Lz_1 = consts_1[1], Q_1 = consts_1[2];
      bool success = helper::kerr_r_limits(Eps_1, Lz_1, Q_1, a_start, M_start, rmin, rmax);
      assert(success);
      success = helper::kerr_th_limits(Eps_1,Lz_1,Q_1, a_start, M_start, thmin, thmax);
      assert(success);
      const double Sr_1 = helper::I_r(Eps_1, Lz_1, Q_1, a_start, M_start, rmin, rmax);
      const double Sth_1= helper::I_th(Eps_1, Q_1, Lz_1, a_start, M_start, thmin, thmax);
      double result = 0;
#ifndef NDEBUG
      Solver solver = Solver::hybrid;
      int growthSuccessful2 = helper::grow_black_hole(consts, a_start, a_2, M_start, M_2, solver, consts_2, result);
      if( growthSuccessful2 == SUCCESSFUL ) traditionalGrowthSuccesses += 1;
#endif
      growthSuccessful = lookup.GrowBlackHole(Sr_1, Sth_1, Lz_1, M_start, a_start, M_2, a_2, consts_2);
#ifndef NDEBUG
      if( growthSuccessful == SUCCESSFUL ) lookupGrowthSuccesses += 1;
      assert( growthSuccessful != SUCCESSFUL && growthSuccessful2 == SUCCESSFUL ); // Lookup table should always perform better, at least in theory
#endif
      
      //cout << "Exiting lookup with growth success, Eps, Lz, Q: " << growthSuccessful << " " << Eps_1 << " " << Lz_1 << " " << Q_1 << " " << consts_2[0] << " " << consts_2[1] << " " << consts_2[2] << endl;
      // Check if the BH growth had a solution:
      if( growthSuccessful == SUCCESSFUL ) {
        // Get orbit:
        {
          const double Eps = consts_2[0];
          const double L = consts_2[1];
          const double Q = consts_2[2];
          if( Eps > 1 || helper::kerr_r_limits(Eps, L, Q, a_2, M_2, rmin,rmax) == false || helper::kerr_th_limits(Eps,L,Q,a_2,M_2,thmin,thmax) == false ) {
            assert(1==1); // This should never happen
            continue;
          }
        }
        // Make sure the particle is sampled on the grid
        if( rmin/M_2 > rMaxLimit ) {
          continue;
        }
        // Particle sampling successful!
        break;
      } else { 
        cout << "Growth unsuccessful :( " << endl;
      }
      if( sampleUnsuccessfulGrowths ) break; // Break anyway
    }
    // Make sure this is a bound orbit:
    const double Eps = consts_2[0];
    const double L = consts_2[1];
    const double Q = consts_2[2];
    assert(helper::kerr_r_limits(Eps, L, Q, a_2, M_2, rmin,rmax) == true && helper::kerr_th_limits(Eps, L, Q, a_2, M_2, thmin,thmax));
    assert(helper::kerr_r_limits(consts_1[0], consts_1[1], consts_1[2], a_start, M_start, rmin,rmax) == true);
    // Sample the particle
    if( growthSuccessful != SUCCESSFUL ) {
      for( int i = 0; i < 4; ++i ) consts_2[i] = numeric_limits<double>::quiet_NaN();
    }
    EpsVectorNew[i] = consts_2[0];
    LzVectorNew[i]  = consts_2[1];
    QVectorNew[i]   = consts_2[2];
    EpsVector[i]    = consts_1[0];
    LzVector[i]     = consts_1[1];
    QVector[i]      = consts_1[2];
    MStartVector[i] = M_start;
    aStartVector[i] = a_start;
    Nin             = Nin+1;
  }
  cout << "Traditional successful growths vs lookup success: " << traditionalGrowthSuccesses << " " << lookupGrowthSuccesses << endl;
}