#ifndef _GRID_H
#define _GRID_H
#include <iostream>
#include <unordered_map>
#include <vector>
#include <array>
#include <cassert>
#include <cmath>
#include <limits>
#include <Eigen/Dense>

# define PI 3.141592653589793238462643383279502884L


typedef double Real;

const Real outOfBounds = std::numeric_limits<double>::quiet_NaN();

using namespace std;


class Grid3D {
public:
  // Initializes grid
  Grid3D( const int radialBins, const int thetaBins, const int phBins, const Real rMax_, const Real rMin_ ) : 
  rDimensions( radialBins ), thDimensions( thetaBins ), phDimensions(phBins), rMax(rMax_), rMin(rMin_), thMax(PI), thMin(0), phMax(2.*PI), phMin(0),  dr((rMax_-rMin_)/(double)(radialBins)), dth((PI-0.)/(double)(thetaBins)), dph((2.*PI-0.)/(double)(phBins)) {
    assert( radialBins>0 && thetaBins>0  && phBins>0 && rMax_> rMin_ && thMax>thMin && phMax>phMin );
    assert( thMin == 0 && thMax == (Real)(PI) );
    assert( phMin == 0 && phMax == (Real)(2.*PI) );
    const int numberOfSpatialGridPoints = rDimensions * thDimensions * phDimensions;
    values.resize( numberOfSpatialGridPoints );
    density.resize( numberOfSpatialGridPoints );
    numberOfSamples.resize( numberOfSpatialGridPoints );
    n.resize( numberOfSpatialGridPoints );
    dT.resize( numberOfSpatialGridPoints );
    numberOfParticles.resize( numberOfSpatialGridPoints );
    for( int i = 0; i < values.size(); ++i ) {values[i] = 0; n[i] = 0; dT[i] = 0; density[i] = 0; numberOfSamples[i] = 0; }
    // Create velocity values:
    velocityValues.reserve( numberOfSpatialGridPoints );
    localVelocityValues.reserve( numberOfSpatialGridPoints);
    localVelocityValues2.reserve(numberOfSpatialGridPoints);
    localMomentumValues.reserve(numberOfSpatialGridPoints);
    localMomentumValues2.reserve(numberOfSpatialGridPoints);
    for( int i = 0; i < values.size(); ++i ) { 
      Eigen::Vector3d tmp(0,0,0); velocityValues.push_back(tmp);
      localVelocityValues.push_back(tmp);
      localVelocityValues2.push_back(tmp);
      localMomentumValues.push_back(tmp);
      localMomentumValues2.push_back(tmp);
      localUnnormalizedDispersions2.push_back(tmp);
      localSecondMoments3.push_back(tmp);
    }
    count = 1;
  }
  
  /*
   * Destructor
   */
  ~Grid3D() {}
  
  /**
   * Does all the sampling on the grid
   * \param t Time coord
   * \param r r coordinate
   * \param th theta coordinate
   * \param ph Phi value
   * \param weight Weight associated with the particle (usually time spent in that grid point)
   **/
  int sample(const Real t, const Real r, const Real th,  Real ph, 
	     const Real vt,const Real vr,const Real vth,const Real vph, const Real weight, const double a, const double M);
  
  /**
   * Does all the sampling on the grid
   * \param t Time coord
   * \param r r coordinate
   * \param th theta coordinate
   * \param ph Phi value
   * \param weight Weight associated with the particle (usually time spent in that grid point)
   **/
  int sample(const vector<Real> & t, const vector<Real> & r, const vector<Real> & th,  vector<Real> & ph, 
	     const vector<Real> & vt,const vector<Real> & vr,const vector<Real> & vth,const vector<Real> & vph, const vector<Real> & weight, const double a, const double M, const double T);

  /**
   * Adds a value to a grid point
   * \param r r coordinate
   * \param th theta coordinate
   * \param value Value to insert
   **/
  void value( const Real r, const Real th, const Real ph, const Real val );
  
  /**
   * Returns the value at a grid point
   * 
   * \param r r coordinate
   * \param th theta coordinate
   **/
  Real value( const Real r, const Real th, const Real ph );

  /**
   * Adds a value to a grid point
   * \param r r coordinate
   * \param th theta coordinate
   * \param value Value to insert
   **/
  void velocityValue( const Real r, const Real th, const Real ph, const Eigen::Vector3d & value, const Real dt );

  /**
   * Adds a value to a grid point
   * \param r r coordinate
   * \param th theta coordinate
   * \param value Value to insert
   **/
  void timeValue( const Real r, const Real th, const Real ph, const double value );

  /**
   * Adds a value to a grid point
   * \param r r coordinate
   * \param th theta coordinate
   * \param value Value to insert
   **/
  double timeValue( const Real r, const Real th, const Real ph );
  
  /**
   * Returns the value at a grid point
   * 
   * \param r r coordinate
   * \param th theta coordinate
   **/
  Eigen::Vector3d velocityValue( const Real r, const Real th, const Real ph );

  /**
   * Adds a value to a grid point
   * \param r r coordinate
   * \param th theta coordinate
   * \param value Value to insert
   **/
  void localVelocityValue( const Real r, const Real th, const Real ph, const Eigen::Vector3d & value, const Real dt );

  /**
   * Returns the value at a grid point
   * 
   * \param r r coordinate
   * \param th theta coordinate
   **/
  Eigen::Vector3d localVelocityValue( const Real r, const Real th, const Real ph );
  
  
  /**
   * Returns the value at a grid point
   * 
   * \param r r coordinate
   * \param th theta coordinate
   **/
  Eigen::Vector3d localVelocityValue2( const Real r, const Real th, const Real ph );


  /**
   * Adds +1 to the counter
   */
  void counter() { count += 1; }

    
  /**
   * Writes out the grid structure
   * 
   * \param fileName Name of the file to write to
   * 
   **/
  void save_grid( const std::string fileName );

  /**
   * Normalizes the grid values with the constant A
   **/
  void normalize_grid( const double NTot );
  
  void normalize_grid( const unsigned int Ntotal, const  unsigned int Nsampled, const double integration_time, const double rho0, const double N0 );

  
  const int rDimensions; const int thDimensions; const int phDimensions; const Real rMax, rMin, thMax, thMin, phMax, phMin;
  const double dr; const double dth; const double dph;

  /*
   * Returns the index in "values" vector corresponding to coordinates e.g. maps (r,theta)->(index), which can be used by values[index]
   **/
  uint32_t return_index( const Real r, const Real th, const Real ph );

  
  /**
   * Returns r and theta from index at the middle of the cell
   **/
  void return_r_theta( const uint32_t index, Real &r, Real &th, Real &ph );
  
  //unordered_map<uint32_t, Real> values;
  vector<Real> values; // 3D real space
  vector<Real> dT; // Total integrated time in a cell
  vector<Real> n; // Number density
  vector<Real> numberOfParticles;
  vector<Eigen::Vector3d> velocityValues; // 3D real space (summed)
  vector<Eigen::Vector3d> localVelocityValues; // velocity values in local frame
  vector<Eigen::Vector3d> localVelocityValues2; // velocity values in local frame
  vector<Eigen::Vector3d> localMomentumValues; // velocity values in local frame
  vector<Eigen::Vector3d> localMomentumValues2; // velocity values in local frame
  vector<Eigen::Vector3d> localUnnormalizedDispersions2; // Dispersions in local frame :)
  vector<Eigen::Vector3d> localSecondMoments3; // Second moment in local frame!
  vector<Real> density;
  vector<Real> numberOfSamples;
  
  unsigned int count;

  ///////
private:
  uint32_t return_r_index( const Real r );
  
  uint32_t return_th_index( const Real th );
  
  uint32_t return_ph_index( const Real ph );
};
#endif
