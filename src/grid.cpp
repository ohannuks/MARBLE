#include "grid.h"
#include <fstream>
#include <cmath>
#include "helperFunctions.h"

using namespace Eigen;

inline 
uint32_t Grid3D::return_r_index(const Real r)
{
  assert( rMin <= r && r <= rMax );
  const int rIndex = (r - rMin) / ((rMax - rMin) / rDimensions);
  assert( rIndex >= 0 );
  return (uint32_t)rIndex;
}

inline 
uint32_t Grid3D::return_th_index(const Real th)
{
  assert( thMin <= th && th <= thMax );
  const int thIndex = (th - thMin) / ((thMax - thMin) / thDimensions);
  assert( thIndex >= 0 );
  return (uint32_t)thIndex;
}

inline 
uint32_t Grid3D::return_ph_index(const Real ph)
{
  assert( phMin <= ph && ph <= phMax );
  const int phIndex = (ph - phMin) / ((phMax - phMin) / phDimensions);
  assert( phIndex >= 0 );
  return (uint32_t)phIndex;
}

inline
uint32_t Grid3D::return_index(const Real r, const Real th, const Real ph)
{
  assert( r <= rMax && th <= thMax && ph <= phMax );
  assert( rMin <= r && thMin <= th && phMin <= ph );
  // If you dont get what is happening, draw a picture ;)
  const uint32_t rIndex  = return_r_index(r);
  const uint32_t thIndex = return_th_index(th);
  const uint32_t phIndex = return_ph_index(ph);
  return rIndex * thDimensions * phDimensions + thIndex * phDimensions + phIndex;
}

void Grid3D::value(const Real r, const Real th, const Real ph_, const Real val)
{
  assert( val >= 0 );
  Real ph = ph_;
  if( ph < phMin || phMax < ph ) { ph = fmod(ph,(2.*M_PI)); if(ph < 0 ) ph = ph + 2.*M_PI; }
  if( th > thMax || th<thMin ) { cout << "Th: " << th << endl; return;}
  if( ph > phMax || ph < phMin ) { cout << "Ph: " << ph << endl; return; }
  if( r < rMin || rMax < r ) { return;} // Out of bounds
  const uint32_t index = return_index(r, th, ph);
  assert( index < values.size() );
  values[index] += val;
  return;
}

void Grid3D::velocityValue(const Real r, const Real th, const Real ph_, const Vector3d & val, const Real dt)
{
  Real ph = ph_;
  if( ph < phMin || phMax < ph ) { ph = fmod(ph,(2.*M_PI)); if(ph < 0 ) ph = ph + 2.*M_PI; }
  if( th > thMax || th<thMin ) { cout << "Th: " << th << endl; return;}
  if( ph > phMax || ph < phMin ) { cout << "Ph: " << ph << endl; return; }
  if( r < rMin || rMax < r ) { return;} // Out of bounds
  const uint32_t index = return_index(r, th, ph);
  assert( index < velocityValues.size() );
  velocityValues[index] += dt*val;
  return;
}

// Input: Velocities
void Grid3D::localVelocityValue(const Real r, const Real th, const Real ph_, const Vector3d & val, const Real dt)
{
  Real ph = ph_;
  if( ph < phMin || phMax < ph ) { ph = fmod(ph,(2.*M_PI)); if(ph < 0 ) ph = ph + 2.*M_PI; }
  if( th > thMax || th<thMin ) { cout << "Th: " << th << endl; return;}
  if( ph > phMax || ph < phMin ) { cout << "Ph: " << ph << endl; return; }
  if( r < rMin || rMax < r ) { return;} // Out of bounds
  const uint32_t index = return_index(r, th, ph);
  assert( index < velocityValues.size() );
  // Add to local velocity:
  localVelocityValues[index] += dt * val;
  localVelocityValues2[index] += dt * (val.array()*val.array()).matrix();
  return;
}

// Input: Velocities
void Grid3D::timeValue(const Real r, const Real th, const Real ph_, const double val)
{
  Real ph = ph_;
  if( ph < phMin || phMax < ph ) { ph = fmod(ph,(2.*M_PI)); if(ph < 0 ) ph = ph + 2.*M_PI; }
  if( th > thMax || th<thMin ) { cout << "Th: " << th << endl; return;}
  if( ph > phMax || ph < phMin ) { cout << "Ph: " << ph << endl; return; }
  if( r < rMin || rMax < r ) { return;} // Out of bounds
  const uint32_t index = return_index(r, th, ph);
  assert( index < dT.size() );
  // Add to local value:
  dT[index] += val;
  return;
}


double Grid3D::timeValue(const Real r, const Real th, const Real ph)
{
  if( th > thMax || th<thMin ) { cout << "Th: " << th << endl; return outOfBounds;}
  if( r < rMin || rMax < r ) return outOfBounds; // Check with isNaN
  const uint32_t index = return_index(r, th,ph);
  const Real value = dT[index];
  assert( value >= 0 );
  return value;
}


Real Grid3D::value(const Real r, const Real th, const Real ph)
{
  if( th > thMax || th<thMin ) { cout << "Th: " << th << endl; return outOfBounds;}
  if( r < rMin || rMax < r ) return outOfBounds; // Check with isNaN
  const uint32_t index = return_index(r, th,ph);
  const Real value = values[index];
  assert( value >= 0 );
  return value;
}


Vector3d Grid3D::velocityValue(const Real r, const Real th, const Real ph)
{
  const Vector3d outOfBoundsVector(outOfBounds,outOfBounds,outOfBounds);
  if( th > thMax || th<thMin ) { cout << "Th: " << th << endl; return outOfBoundsVector;}
  if( r < rMin || rMax < r ) return outOfBoundsVector; // Check with isNaN
  const uint32_t index = return_index(r, th,ph);
  const Vector3d velocityValue = velocityValues[index];
  return velocityValue;
}

Vector3d Grid3D::localVelocityValue( const Real r, const Real th, const Real ph ) {
  const Vector3d outOfBoundsVector(outOfBounds,outOfBounds,outOfBounds);
  if( th > thMax || th<thMin ) { cout << "Th: " << th << endl; return outOfBoundsVector;}
  if( r < rMin || rMax < r ) return outOfBoundsVector; // Check with isNaN
  const uint32_t index = return_index(r, th,ph);
  const Vector3d velocityValue = localVelocityValues[index];
  return velocityValue;
}

Vector3d Grid3D::localVelocityValue2( const Real r, const Real th, const Real ph ) {
  const Vector3d outOfBoundsVector(outOfBounds,outOfBounds,outOfBounds);
  if( th > thMax || th<thMin ) { cout << "Th: " << th << endl; return outOfBoundsVector;}
  if( r < rMin || rMax < r ) return outOfBoundsVector; // Check with isNaN
  const uint32_t index = return_index(r, th,ph);
  const Vector3d velocityValue = localVelocityValues2[index];
  return velocityValue;
}


void Grid3D::return_r_theta(const uint32_t index, Real& r, Real& th, Real& ph)
{
  assert( index < values.size() ); // Cant be bigger than values container
  const uint32_t rIndex = index/(thDimensions*phDimensions); // Get rIndex
  const uint32_t thIndex = (index/phDimensions)%thDimensions;
  const uint32_t phIndex = index%phDimensions; // Get phIndex
  // Input values:
  r = rMin   + (Real)(rMax-rMin)  / (Real)(rDimensions)  * (rIndex +0.5);
  th = thMin + (Real)(thMax-thMin)/ (Real)(thDimensions) * (thIndex+0.5);
  ph = phMin + (Real)(phMax-phMin)/ (Real)(phDimensions) * (phIndex+0.5);
  return;
}


void Grid3D::save_grid(const string fileName)
{
  ofstream file;
  
  file.open(fileName, ios::out | ios::binary);
  
  // Output into binary:
  const int numberOfVariables = 20;
  const int numberOfAttributes = 11;
  float attributes[numberOfAttributes] = {(float)rMin, (float)thMin, (float)phMin, (float)rMax, (float)thMax, (float)phMax,(float)rDimensions, (float)thDimensions, (float)phDimensions, (float)values.size(), (float)numberOfVariables};
  file.write((char*)attributes, numberOfAttributes*sizeof(float)); // Write attributes
  Real r, th, ph;
  //file << rMin << " " << thMin << " " << phMin << " " << 0 << endl;
  //file << rMax << " " << thMax << " " << phMax << " " << 0 << endl;
  //file << rDimensions << " " << thDimensions << " " << phDimensions << " " << 0 << endl;
  vector<float> rVals, thVals, phVals, vals, vrVals,vthVals,vphVals, timeVals, densityVals, numberOfSamplesVals, numberOfParticlesVals;
  rVals.resize(values.size());
  thVals.resize(values.size());
  phVals.resize(values.size());
  vals.resize(values.size());
  vrVals.resize(values.size());
  vthVals.resize(values.size());
  vphVals.resize(values.size());
  timeVals.resize(values.size());

  // Local velocity values:
  vector<float> Vx, Vy, Vz, Vx_sigma, Vy_sigma, Vz_sigma;
  Vx.resize(values.size());
  Vy.resize(values.size());
  Vz.resize(values.size());
  Vx_sigma.resize(values.size());
  Vy_sigma.resize(values.size());
  Vz_sigma.resize(values.size());
  vector<float> Px, Py, Pz, Px_sigma, Py_sigma, Pz_sigma;
  Px.resize(values.size());
  Py.resize(values.size());
  Pz.resize(values.size());
  Px_sigma.resize(values.size());
  Py_sigma.resize(values.size());
  Pz_sigma.resize(values.size());
  densityVals.resize(values.size());
  numberOfSamplesVals.resize(values.size());
  numberOfParticlesVals.resize(values.size());

  for( int i = 0; i < values.size(); ++i ) {
    return_r_theta((uint32_t)i, r, th,ph); // Get r, theta
    rVals[i] = (float)r;
    thVals[i] = (float)th;
    phVals[i] = (float)ph;
    const uint32_t index = return_index(r, th,ph);
    vals[i] = (float)value(r, th,ph);
    timeVals[i] = (float)timeValue(r,th,ph);
    double totalTime = timeVals[i];
    if( totalTime == 0 ) totalTime = numeric_limits<double>::max();
    // Get velocity value:
    const Vector3d velocityVals = velocityValue(r,th,ph);
    //file << r << " " << th << " " << ph << " " << value(r, th,ph) << endl;
    vrVals[i]  = velocityVals[0]/(double)(totalTime);
    vthVals[i] = velocityVals[1]/(double)(totalTime);
    vphVals[i] = velocityVals[2]/(double)(totalTime);
    // Get local velocity values:
    const Vector3d localVelocityVals = localVelocityValue(r,th,ph);
    const Vector3d localVelocityVals2 = localVelocityValue2(r,th,ph);
    Vx[i] = localVelocityVals[0] / (double)(totalTime);
    Vy[i] = localVelocityVals[1] / (double)(totalTime);
    Vz[i] = localVelocityVals[2] / (double)(totalTime);
    Vx_sigma[i] = localVelocityVals2[0]/(double)(totalTime)-pow(localVelocityVals[0]/(double)(totalTime),2);
    Vy_sigma[i] = localVelocityVals2[1]/(double)(totalTime)-pow(localVelocityVals[1]/(double)(totalTime),2);
    Vz_sigma[i] = localVelocityVals2[2]/(double)(totalTime)-pow(localVelocityVals[2]/(double)(totalTime),2);
    if( !((Vx_sigma[i] >= 0 && Vy_sigma[i] >= 0 && Vz_sigma[i] >= 0)) ) cout << Vx_sigma[i] << " " << Vy_sigma[i] << " " << Vz_sigma[i] << endl;
    //assert(Vx_sigma[i] >= 0 && Vy_sigma[i] >= 0 && Vz_sigma[i] >= 0);
    // Save also local momentum:
    const Vector3d localMomentumVals = localMomentumValues[index];
    const Vector3d localMomentumVals2 = localMomentumValues2[index];
    Px[i] = localMomentumVals[0] / (double)(totalTime);
    Py[i] = localMomentumVals[1] / (double)(totalTime);
    Pz[i] = localMomentumVals[2] / (double)(totalTime);
    Px_sigma[i] = localMomentumVals2[0]/(double)(totalTime)-pow(localMomentumVals[0]/(double)(totalTime),2);
    Py_sigma[i] = localMomentumVals2[1]/(double)(totalTime)-pow(localMomentumVals[1]/(double)(totalTime),2);
    Pz_sigma[i] = localMomentumVals2[2]/(double)(totalTime)-pow(localMomentumVals[2]/(double)(totalTime),2);
    densityVals[i] = density[index];
    numberOfSamplesVals[i] = numberOfSamples[index];
    numberOfParticlesVals[i] = numberOfParticles[index];
  }
  file.write((char*)&rVals[0], rVals.size()*sizeof(float));
  file.write((char*)&thVals[0], thVals.size()*sizeof(float));
  file.write((char*)&phVals[0], phVals.size()*sizeof(float));
  file.write((char*)&vals[0], vals.size()*sizeof(float));
  file.write((char*)&vrVals[0],  vals.size()*sizeof(float));
  file.write((char*)&vthVals[0], vals.size()*sizeof(float));
  file.write((char*)&vphVals[0], vals.size()*sizeof(float));
  file.write((char*)&Vx[0], vals.size()*sizeof(float));
  file.write((char*)&Vy[0], vals.size()*sizeof(float));
  file.write((char*)&Vz[0], vals.size()*sizeof(float));
  file.write((char*)&Vx_sigma[0], vals.size()*sizeof(float));
  file.write((char*)&Vy_sigma[0], vals.size()*sizeof(float));
  file.write((char*)&Vz_sigma[0], vals.size()*sizeof(float));
  file.write((char*)&timeVals[0], vals.size()*sizeof(float));
  file.write((char*)&Px[0], vals.size()*sizeof(float));
  file.write((char*)&Py[0], vals.size()*sizeof(float));
  file.write((char*)&Pz[0], vals.size()*sizeof(float));
  file.write((char*)&Px_sigma[0], vals.size()*sizeof(float));
  file.write((char*)&Py_sigma[0], vals.size()*sizeof(float));
  file.write((char*)&Pz_sigma[0], vals.size()*sizeof(float));
  file.write((char*)&densityVals[0], vals.size()*sizeof(float));
  file.write((char*)&numberOfSamplesVals[0], vals.size()*sizeof(float));
  file.write((char*)&numberOfParticlesVals[0], vals.size()*sizeof(float));
  file.close();
  return;
}

// Record *everything*
// OUTDATED
int Grid3D::sample(const Real t, const Real r, const Real th,  Real ph, 
                   const Real vt,const Real vr,const Real vth,const Real vph, const Real weight, const double a, const double M) {
  if( ph < phMin || phMax < ph ) { ph = fmod(ph,(2.*M_PI)); if(ph < 0 ) ph = ph + 2.*M_PI; }
  if( r > rMax || th > thMax || ph > phMax ) return UNSUCCESSFUL;
  if( r < rMin || th < thMin || ph < phMin ) return UNSUCCESSFUL;
  //4-pos, vel
  const double x[4] = {t,r,th,ph};
  const double u[4] = {vt,vr,vth,vph};
  // Get metric:
  const double a2 = a*a;
  const double M2 = M*M;
  const double costh = cos(th); const double costh2 = costh*costh;
  const double sinth = sin(th); const double sinth2 = sinth*sinth;
  const double r2 = r*r;
  double g[4][4];
  helper::gmunu(x,a,M,g);
  // Get determinant of metric:
  // NOTE: This is incorrect!FIX!
  const double detg = -1.*pow(a2*costh2+r2,2)*sinth2;
  // Get index:
  const int index = return_index(r,th,ph);
  // Record *all* the quantities:
  // Density in coordinate frame:
  const double density = weight / sqrt(-1.*detg);
  // Velocity
  const Eigen::Vector3d v(vr/vt,vth/vt,vph/vt);
  // Transform velocity to local frame:
  double tetrad[4][4] = {0};
  helper::ZAMOTetrad(x, a, M, tetrad, true);
  double localVel[4] = {0};
  for( int m = 0; m < 4; ++m ) for( int mu = 0; mu < 4; ++mu ) {
    localVel[m] += tetrad[m][mu]*u[mu];
  }
  // Make sure its physical:
  const double c = 1;
  for(int i = 1; i < 4; ++i ) { if( localVel[i]/localVel[0] > c) return UNSUCCESSFUL; }
  /// Calculate velocity in local tetrad frame:
  const Eigen::Vector3d localv(localVel[1]/localVel[0], localVel[1]/localVel[0], localVel[2]/localVel[0]);
  // Calculate velocity pow 2
  const Eigen::Vector3d localv2(localv(0)*localv(0),localv(1)*localv(1),localv(2)*localv(2));
  /// Calculate velocity in local tetrad frame:
  const Eigen::Vector3d localMomentum(localVel[1], localVel[1], localVel[2]);
  // Calculate velocity pow 2
  const Eigen::Vector3d localMomentum2(localMomentum(0)*localMomentum(0),localMomentum(1)*localMomentum(1),localMomentum(2)*localMomentum(2));

  n[index] += density;
  values[index] += density;
  dT[index] += weight;
  cout << "I am adding " << weight << " to numberOfParticles " << numberOfParticles[index] << endl;
  numberOfParticles[index] += weight;
  velocityValues[index] += weight*v;
  localVelocityValues[index] += weight*localv;
  localVelocityValues2[index]+= weight*localv2;
  localMomentumValues[index] += weight*localMomentum;
  localMomentumValues2[index]+= weight*localMomentum2;
  return SUCCESSFUL;
}

void Grid3D::normalize_grid(const double NTot)
{
  for( int i = 0 ; i < values.size(); ++i ) {
    density[i] = values[i] / NTot;
    numberOfParticles[i] = numberOfParticles[i] / NTot;
  }
//  const double pi = M_PI;
//  const double c = 63198; // AU/yr
//  const double G = 4.*pi*pi; // AU^3/(yr^2 Msun)
//  const double RPhysical = MPhysical*G/(c*c);
//  const double toGeVPerCmCubic = 3.331702e17;
//  for( int i = 0; i < values.size(); ++i ) {
//    density[i] = values[i] * A / pow(RPhysical,3); // Msun / AU^3
//    density[i] = values[i]*toGeVPerCmCubic; // GeV/cm^3
//  }
  return;
}

void Grid3D::normalize_grid( const unsigned int Ntotal, const unsigned int Nsampled, const double integration_time, const double rho0, const double N0 ) {
  // Solve Nenclosed:
  const double Nenc = Nsampled;
  // Solve Menclosed
  const double Menc = ( (double)Nenc / (double)Ntotal ) * rho0/N0 * integration_time;
}

// Record *everything*
int Grid3D::sample(const vector<Real> & t_, const vector<Real> & r_, const vector<Real> & th_,  vector<Real> & ph_, 
           const vector<Real> & vt_,const vector<Real> & vr_,const vector<Real> & vth_,const vector<Real> & vph_, const vector<Real> & weight_, 
           const double a, const double M, const double T) {
  assert( t_.size() == vr_.size() );
  vector<Real> samples; samples.resize(values.size()); for( int i = 0; i < values.size(); ++i ) samples[i] = 0;
  const int nElements = t_.size();
  // Create tetrad frame vectors:
  vector<Real> ut_, ur_, uth_, uph_;
  ut_.resize(nElements); ur_.resize(nElements); uth_.resize(nElements); uph_.resize(nElements);
  for( int i = 0; i < nElements; ++i ) {
    const double t  = t_[i];
    const double r  = r_[i];
    const double th = th_[i];
    double ph = ph_[i];
    const double vt  = vt_[i];
    const double vr  = vr_[i];
    const double vth = vth_[i];
    const double vph = vph_[i];
    const double x[4] = {t,r,th,ph};
    const double v[4] = {vt,vr,vth,vph};

    if( ph < phMin || phMax < ph ) { ph = fmod(ph,(2.*M_PI)); if(ph < 0 ) ph = ph + 2.*M_PI; }
    if( r > rMax || th > thMax || ph > phMax ) return UNSUCCESSFUL;
    if( r < rMin || th < thMin || ph < phMin ) return UNSUCCESSFUL;
    // Make sure the velocity is physical
    // 4-pos, vel
    {
      double g[4][4] = {0};
      helper::gmunu(x,a,M,g);
      double normalization = 0;
      for( int mu = 0; mu < 4; ++mu ) for( int nu = 0; nu < 4; ++nu ) { normalization += g[mu][nu]*v[mu]*v[nu]; } // Should be -1
      if( abs( normalization + 1 ) > 0.0001 ) return UNSUCCESSFUL;
    }

    // Transform velocity to local frame:
    double tetrad[4][4] = {0};
    helper::ZAMOTetrad(x, a, M, tetrad, true);
    double u[4] = {0};
    for( int m = 0; m < 4; ++m ) for( int mu = 0; mu < 4; ++mu ) {
      u[m] += tetrad[m][mu]*v[mu];
    }
    // Make sure the velocity is physical:
    {
      double normalization = -1.*u[0]*u[0] + u[1]*u[1] + u[2]*u[2] + u[3]*u[3];
      if( abs( normalization + 1 ) > 0.0001 ) return UNSUCCESSFUL;
    }
    // Record the velocity:
    ut_[i] = u[0];
    ur_[i] = u[1];
    uth_[i]= u[2];
    uph_[i]= u[3];
  }
  // Save to grid:
  for( int i = 0; i < nElements; ++i ) {
    const double t  = t_[i];
    const double r  = r_[i];
    const double th = th_[i];
    double ph = ph_[i];
    if( ph < phMin || phMax < ph ) { ph = fmod(ph,(2.*M_PI)); if(ph < 0 ) ph = ph + 2.*M_PI; }
    const double vt  = vt_[i];
    const double vr  = vr_[i];
    const double vth = vth_[i];
    const double vph = vph_[i];
    const double weight= weight_[i] / T; // Equals w*dt / T where w is traditionally one, but under adiabatic growth might include thedegeneracy factor
    const double x[4] = {t,r,th,ph};
    // Get metric:
    const double a2 = a*a;
    const double M2 = M*M;
    const double costh = cos(th); const double costh2 = costh*costh;
    const double sinth = sin(th); const double sinth2 = sinth*sinth;
    const double r2 = r*r;
    double g[4][4];
    helper::gmunu(x,a,M,g);
    // Get determinant of metric:
    const double detg = -1.*pow(a2*costh2+r2,2)*sinth2;
    // Get index:
    const int index = return_index(r,th,ph);
    // Record *all* the quantities:
    // Density in coordinate frame:
    const double density = weight / (sqrt(-1.*detg)*dr*dth*dph); // TODO: Should there be dr dth dphi here? Also, what unit is this exactly?
    // Velocity
    const Eigen::Vector3d v(vr/vt,vth/vt,vph/vt);
    // Get the local velocity frame:
    const double localVel[4] = {ut_[i], ur_[i], uth_[i], uph_[i]};
    /// Calculate velocity in local tetrad frame:
    const Eigen::Vector3d localv(localVel[1]/localVel[0], localVel[2]/localVel[0], localVel[3]/localVel[0]);
    // Calculate velocity pow 2
    const Eigen::Vector3d localv2(localv(0)*localv(0),localv(1)*localv(1),localv(2)*localv(2));
    /// Calculate velocity in local tetrad frame:
    const Eigen::Vector3d localMomentum(localVel[1], localVel[2], localVel[3]);
    // Calculate velocity pow 2
    const Eigen::Vector3d localMomentum2(localMomentum(0)*localMomentum(0),localMomentum(1)*localMomentum(1),localMomentum(2)*localMomentum(2));
    
    n[index] += density;
    values[index] += density;
    dT[index] += weight;
    numberOfParticles[index] += weight;
    velocityValues[index] += weight*v;
    localVelocityValues[index] += weight*localv;
    localVelocityValues2[index]+= weight*localv2;
    localMomentumValues[index] += weight*localMomentum;
    localMomentumValues2[index]+= weight*localMomentum2;
    samples[index] = 1;
  }
  for( int i = 0; i < samples.size(); ++i ) {
    numberOfSamples[i] += samples[i];
  }
  return SUCCESSFUL;
}















