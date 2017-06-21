
#include "sampling.h"
#include "Eigen/Dense"
#include "cubature/cubature.h"

const double pi = M_PI;

using namespace std;
using namespace Eigen;

namespace sampling {

int maxwellian_sampling( normal_distribution<double> & normal_dist, mt19937_64 & gen, const double mean_v[3], const double x[4], double v[3] ) {
  while( true ) {
    const double velocity_variance = normal_dist( gen ); // Get the variance
    const double velocity_variance_r = normal_dist( gen ); // Get the variance
    const double velocity_variance_th = normal_dist( gen ); // Get the variance
    const double velocity_variance_ph = normal_dist( gen ); // Get the variance
    if( velocity_variance >= 1 ) continue;
    if( velocity_variance_r >= 1. || velocity_variance_ph >= 1. || velocity_variance_th >= 1. ) continue;
    // theta and phi
    uniform_real_distribution<> costheta(-1.,1.);
    uniform_real_distribution<> phi(0,2*pi);
    //const double th = theta(gen);
    const double ph = phi(gen);
    const double th = acos(costheta(gen));
    // Input vs in cartesian
    const double vv = velocity_variance;
    v[0] = mean_v[0] + vv*sin(th)*cos(ph); v[0] = mean_v[0]+velocity_variance_r;
    v[1] = mean_v[1] + vv*sin(th)*sin(ph); v[1] = mean_v[1]+velocity_variance_th;
    v[2] = mean_v[2] + vv*cos(th); v[2] = mean_v[2]+velocity_variance_ph;
    
    if( v[0]*v[0]+v[1]*v[1]+v[2]*v[2] < 0.95 ) break;
  }
  // Transform to {dr/dt, d\theta/dt, d\phi/dt}
  v[0] = v[0];
  v[1] = v[1]/x[1]; // v=r d\theta/dt
  v[2] = v[2]/(x[1]*sin(x[2])); // v=r d\phi / dt
  return SUCCESSFUL;
}
  
// Sample the r_influence distribution:
int r_sampling( const double rMinInfluence, const double rMaxInfluence, mt19937_64 & gen, double x[4], Profile profile, const double rMaxSample ) {
  // random uniform sampling:
  double r = 0;
  double A = 0;
  const int n = profile;
  assert( Profile::r2 == 2);
  if( profile == Profile::r3 ) {
    A = 1./log(rMaxInfluence/rMinInfluence);
  } else {
    A = ((-3. + n)*pow(rMaxInfluence*rMinInfluence,n))/(pow(rMaxInfluence,n)*pow(rMinInfluence,3) - pow(rMaxInfluence,3)*pow(rMinInfluence,n));
  }
  // rMaxSample simply gives us the ability to optimize u
  uniform_real_distribution<> u(0.,(A*(pow(rMaxSample,3 - n) - pow(rMinInfluence,3 - n)))/(3 - n));
  r = pow((A*pow(rMinInfluence,3 - n) + (3. - n)*u(gen))/A,1./(3. - n));
  // Randomly sample also th and ph
  uniform_real_distribution<> costheta(-1.,1.);
  uniform_real_distribution<> phi(0,2*pi);
  const double th = acos(costheta(gen));
  const double ph = phi(gen);
  x[1] = r;
  x[2] = th;
  x[3] = ph;
  return SUCCESSFUL;
}


int f(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) {
  double sigma = ((double *) fdata)[0];
  double rInfluence = ((double *) fdata)[1];
  double a = ((double*)fdata)[2];
  double M = ((double*)fdata)[3];
  double sum = 0;
  unsigned i;
  for (i = 0; i < ndim; ++i) sum += x[i] * x[i];
  // compute the output value: note that fdim should == 1 from below
  fval[0] = (4.*pi*sum)/pow(sqrt(2.*pi*sigma*sigma),3) * exp(-.5*sum/(sigma*sigma));
  return 0; // success
}











}


