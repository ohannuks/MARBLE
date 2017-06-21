#include <cmath>
//#include "libcmaes/cmaes.h"
#include "rpoly_ak1.h"
#include <algorithm>
#include <gsl/gsl_integration.h>
#include <cassert>
#include <iostream>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>
#include <limits>
#include <Eigen/Dense>
#include "helperFunctions.h"
#include "cubature/cubature.h"
#include <omp.h>
#include <random>
#include <gsl/gsl_siman.h>


#define SUCCESSFULTOLERANCE 1.e-6
#define FLOAT_ERR 1.e-6
#define RHARDCUTOFF 800.0


// Matrix inversion :)
bool gluInvertMatrix(const double m[16], double invOut[16])
{
    double inv[16], det;
    int i;

    inv[0] = m[5]  * m[10] * m[15] - 
             m[5]  * m[11] * m[14] - 
             m[9]  * m[6]  * m[15] + 
             m[9]  * m[7]  * m[14] +
             m[13] * m[6]  * m[11] - 
             m[13] * m[7]  * m[10];

    inv[4] = -m[4]  * m[10] * m[15] + 
              m[4]  * m[11] * m[14] + 
              m[8]  * m[6]  * m[15] - 
              m[8]  * m[7]  * m[14] - 
              m[12] * m[6]  * m[11] + 
              m[12] * m[7]  * m[10];

    inv[8] = m[4]  * m[9] * m[15] - 
             m[4]  * m[11] * m[13] - 
             m[8]  * m[5] * m[15] + 
             m[8]  * m[7] * m[13] + 
             m[12] * m[5] * m[11] - 
             m[12] * m[7] * m[9];

    inv[12] = -m[4]  * m[9] * m[14] + 
               m[4]  * m[10] * m[13] +
               m[8]  * m[5] * m[14] - 
               m[8]  * m[6] * m[13] - 
               m[12] * m[5] * m[10] + 
               m[12] * m[6] * m[9];

    inv[1] = -m[1]  * m[10] * m[15] + 
              m[1]  * m[11] * m[14] + 
              m[9]  * m[2] * m[15] - 
              m[9]  * m[3] * m[14] - 
              m[13] * m[2] * m[11] + 
              m[13] * m[3] * m[10];

    inv[5] = m[0]  * m[10] * m[15] - 
             m[0]  * m[11] * m[14] - 
             m[8]  * m[2] * m[15] + 
             m[8]  * m[3] * m[14] + 
             m[12] * m[2] * m[11] - 
             m[12] * m[3] * m[10];

    inv[9] = -m[0]  * m[9] * m[15] + 
              m[0]  * m[11] * m[13] + 
              m[8]  * m[1] * m[15] - 
              m[8]  * m[3] * m[13] - 
              m[12] * m[1] * m[11] + 
              m[12] * m[3] * m[9];

    inv[13] = m[0]  * m[9] * m[14] - 
              m[0]  * m[10] * m[13] - 
              m[8]  * m[1] * m[14] + 
              m[8]  * m[2] * m[13] + 
              m[12] * m[1] * m[10] - 
              m[12] * m[2] * m[9];

    inv[2] = m[1]  * m[6] * m[15] - 
             m[1]  * m[7] * m[14] - 
             m[5]  * m[2] * m[15] + 
             m[5]  * m[3] * m[14] + 
             m[13] * m[2] * m[7] - 
             m[13] * m[3] * m[6];

    inv[6] = -m[0]  * m[6] * m[15] + 
              m[0]  * m[7] * m[14] + 
              m[4]  * m[2] * m[15] - 
              m[4]  * m[3] * m[14] - 
              m[12] * m[2] * m[7] + 
              m[12] * m[3] * m[6];

    inv[10] = m[0]  * m[5] * m[15] - 
              m[0]  * m[7] * m[13] - 
              m[4]  * m[1] * m[15] + 
              m[4]  * m[3] * m[13] + 
              m[12] * m[1] * m[7] - 
              m[12] * m[3] * m[5];

    inv[14] = -m[0]  * m[5] * m[14] + 
               m[0]  * m[6] * m[13] + 
               m[4]  * m[1] * m[14] - 
               m[4]  * m[2] * m[13] - 
               m[12] * m[1] * m[6] + 
               m[12] * m[2] * m[5];

    inv[3] = -m[1] * m[6] * m[11] + 
              m[1] * m[7] * m[10] + 
              m[5] * m[2] * m[11] - 
              m[5] * m[3] * m[10] - 
              m[9] * m[2] * m[7] + 
              m[9] * m[3] * m[6];

    inv[7] = m[0] * m[6] * m[11] - 
             m[0] * m[7] * m[10] - 
             m[4] * m[2] * m[11] + 
             m[4] * m[3] * m[10] + 
             m[8] * m[2] * m[7] - 
             m[8] * m[3] * m[6];

    inv[11] = -m[0] * m[5] * m[11] + 
               m[0] * m[7] * m[9] + 
               m[4] * m[1] * m[11] - 
               m[4] * m[3] * m[9] - 
               m[8] * m[1] * m[7] + 
               m[8] * m[3] * m[5];

    inv[15] = m[0] * m[5] * m[10] - 
              m[0] * m[6] * m[9] - 
              m[4] * m[1] * m[10] + 
              m[4] * m[2] * m[9] + 
              m[8] * m[1] * m[6] - 
              m[8] * m[2] * m[5];

    det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

    if (det == 0)
        return false;

    det = 1.0 / det;

    for (i = 0; i < 16; i++)
        invOut[i] = inv[i] * det;

    return true;
}

const double pi = M_PI;
const double G = 1;
const int iterations=100000;

using namespace std;

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

namespace helper {
// Find u^0 from u^\mu u_\mu = -1 for converting between newtonian and non-newtonian approach
void findVt( const double x[4], const double v[3], const double a, const double M, double u[4] ) {

  u[1] = v[0];
  u[2] = v[1];
  u[3] = v[2];
  // Find metric
  double g[4][4] = {0};
  helper::gmunu( x,
                 a,
                 M,
                 g );
  long double norm3 = 0;
  for( int mu = 1; mu < 4; ++mu ) for( int nu = 1; nu < 4; ++nu ) norm3 += g[mu][nu] * u[mu]*u[nu];
  // Calculate u^0 from g_00(u^0)^2 + 2 g_03 u^3 u^0 + norm3 = -1
  const long double a_ = g[0][0];
  const long double b = 2.*g[0][3]*u[3];
  const long double c = norm3 + 1;
  u[0] = (float)((-1.*b - sqrt(b*b-4.*a_*c))/(2.*a_));
  assert( u[0] > 0 );
  return;
}

// e_m^\mu
// Transforms velocity into zamo tetrad
void ZAMOTetrad( const double x[4],
            const double a,
            const double M, 
	    double tetrad[4][4], 
	    bool inverse ) {
  const long double r = x[1];
  const long double r2 = r*r;
  const long double th = x[2];
  // Get sin^2 th and cos^2 th
  const long double sinth = sin(th);
  const long double costh = cos(th);
  const long double sinth2=sinth*sinth; 
  const long double costh2=costh*costh;
  const long double a2 = a*a;
  const long double Sigma = r2 + a2 * costh2;
  const long double Delta = r2 - 2.*M*r + a2;
  const long double A = pow(r2+a2,2)-Delta*a2*sinth2;
  const long double gtt=-1.+(2.*M*r)/(a2*costh2+r2);
  const long double gpp=(sinth2*(pow(a2+r2,2)-a2*(a2+r*(-2.*M+r))*sinth2))/(a2*costh2+r2);
  const long double gtp=-1.*((2.*a*M*r*sinth2)/(a2*costh2+r2));
  const long double Omega=-gtp/gpp;
  const long double Omega2=Omega*Omega;

  // Zero out
  for (int i = 0; i < 4; ++i ) for( int j = 0; j < 4; ++j ) {tetrad[i][j]=0;}

  //tetrad[0][0] = (r2+a2)/sqrt(Delta*Sigma);
  //tetrad[1][1] = sqrt(Delta/Sigma);
  //tetrad[2][2] = 1./sqrt(Sigma);
  //tetrad[3][3] = 1./(sqrt(Sigma)*sinth);
  //tetrad[0][3] = a / sqrt(Delta*Sigma);
  //tetrad[3][0] = a * sinth / sqrt(Sigma);
  tetrad[0][0] = 1./sqrt(abs(gtt-Omega2*gpp));
  tetrad[0][3] = Omega/sqrt(abs(gtt-Omega2*gpp));
  tetrad[1][1] = sqrt(Delta/Sigma);
  tetrad[2][2] = 1./sqrt(Sigma);
  tetrad[3][3] = 1./sqrt(gpp);

  if( inverse ) {
    double inverse[4][4] = {0};
    gluInvertMatrix(&tetrad[0][0],&inverse[0][0]);
    for( int i = 0; i < 4; ++i ) for( int j = 0; j < 4; ++j ) {
      tetrad[i][j] = inverse[j][i];
    }
  }
#ifdef TEST_RUN
  if( !inverse ) {
    double eta[4][4] = {0};
    double g[4][4]   = {0};
    gmunu( x, a, M, g );
    for( int m=0; m<4; ++m ) for( int n=0; n<4; ++n ) for( int mu = 0; mu < 4; ++mu ) for( int nu = 0; nu < 4; ++nu ) {
      eta[m][n]=tetrad[m][mu]*tetrad[n][nu]*g[mu][nu];
    }
    // Check off-diagonal and on-diagonal components:
    double offdiagonal = 0;
    double ondiagonal  = 0;
    for( int m = 0; m < 4; ++m ) for( int n = 0; n < 4; ++n ) {
      if( m == n ) ondiagonal += abs(eta[m][n]);
      if( m != n ) offdiagonal+= abs(eta[m][n]);
    }
    if( offdiagonal > 1e-8 ) { cerr << "BAD OFF-DIAGONAL COMPONENT: " << offdiagonal << " (SHOULD BE ZERO) " <<endl; exit(1); }
    if( ondiagonal-4> 1e-8 ) { cerr << "BAD ON-DIAGONAL  COMPONENT: " << ondiagonal  << " (SHOULD BE 4) " <<endl; exit(1); }
  }
#endif
  return;
}

// Makes sure that our solver is being used properly (the solver assumes this condition to ensure working)
bool check_conditions( const double M_1, const double M_2, const double a_1, const double a_2 ) {
  // The ratio of a_1/M_1 and a_2/M_2 is assumed to be the same:
  if( abs(M_1*a_2 - M_2*a_1) < 1.e-10  ) return true;
  return false;
}


void gmunu( const double x[4],
            const double a,
            const double M,
            double g[4][4] ) {
  const double a2=a*a;
  const double t = x[0];
  const double r = x[1];
  const double th= x[2];
  const double ph= x[3];
  const double sth2=pow(sin(th),2);
  const double cth2=pow(cos(th),2);
  const double r2=r*r;
  const double sigma=r2+a2*cth2;
  const double delta=r2-2.*M*r+a2;
  const double gtt=(-(delta/sigma) + (a2*sth2)/sigma);
  const double grr=sigma/delta;
  const double gthth=sigma;
  const double gpp=(sth2*(pow(a2 + r2,2) - a2*delta*sth2))/sigma;
  const double gtp=(-((2.*a*(a2 - delta + r2)*sth2)/sigma))/2.;
  for( int i = 0; i < 4; ++i ) for( int j = 0; j < 4; ++j ) g[i][j]=0;
  g[0][0] = gtt;
  g[1][1] = grr;
  g[2][2] = gthth;
  g[3][3] = gpp;
  g[0][3] = gtp;
  g[3][0] = gtp;
  return;
}

static
double Csc( const double x ) {
  return -2.*sin(x)/(cos(2.*x)-1.);
}

static 
double _R(const double Eps,const double L,const double Q,const double a,const double M,const double r) {
  const double Eps2=Eps*Eps, L2=L*L, a2=a*a;
  const double c4 = (Eps2-1.),                     // *r^4
               c3 = 2.*M,                          // *r^3
               c2 = -1.*L2-Q+a2*(Eps-1.)*(Eps+1.), // *r^2
               c1 = 2.*M*(Q+pow(L-a*Eps,2)),
               c0 = -1.*a2*Q;
  const double r2=r*r, r3=r2*r, r4=r2*r2;
  return c4*r4+c3*r3+c2*r2+c1*r+c0;
}

static 
double _Theta(const double Eps,const double L,const double Q,const double a,const double M,const double th) {
  return pow(a*Eps - L,2) + Q - pow(a,2)*pow(cos(th),2) - pow(-(L*Csc(th)) + a*Eps*sin(th),2);
}

void constants(const double pos[4],
               const double u[4],
               const double a,
               const double M,
               double consts[4]) {
  if( u[0] != u[0] ) {cout << "u0, " << u[0] << "Is nan at " << __FILE__ << " " << __LINE__ << endl;}
  const double a2=a*a;
  double g[4][4] = {0};
  gmunu(pos,a,M,g);
  const double gtt=g[0][0];
  const double grr=g[1][1];
  const double gthth=g[2][2];
  const double gpp=g[3][3];
  const double gtp=g[0][3];
  // Constants of motion:
  const double r = pos[1];
  const double th = pos[2];
  const double sinth2=pow(sin(th),2); const double sinth = sin(th);
  const double costh2=pow(cos(th),2); const double costh = cos(th);
  const double r2 = r*r;
  double Eps = -gtt*u[0]-gtp*u[3]; double Eps2 = Eps*Eps;
  double L = gtp*u[0]+gpp*u[3]; double L2 = L*L;
  //double Q = Omega + costh2*(a2*(1.-Eps2)+L2/sinth2);
  double u_d[4] = {
    gtt*u[0]+gtp*u[3], 
    grr*u[1],
    gthth*u[2],
    gtp*u[0]+gpp*u[3]
  };
  // Eq 8.7.3. claimedBest.pdf
  double division = u_d[3]/sinth;
  if( sinth == 0 && u_d[3] == 0 ) division = 1;
  double K = pow(u_d[0]*a*sinth+u_d[3]/sinth,2)+u_d[2]*u_d[2]-a2*costh2*(-1.);
  // Eq. 8.5.19
  double Q = K-pow(Eps*a-L,2);
  consts[0]=Eps;
  consts[1]=L;
  consts[2]=Q;
  consts[3]=K;
  return;
}


bool kerr_th_limits(const double Eps,
                    const double L,
                    const double Q,
                    const double a,
                    const double M,
                    double & thMin,
                    double & thMax ) {
  const double a2=a*a;

  // Solve using bisection method:
  /////////
  const double tolerance=1.e-14;
  double thLimits[2] = {0};
  for( int i = 0; i < 2; ++i ) {
    double a,b;
    if( i == 0 ) {
      a = 1.e-12;
      b = pi/2.;
    } else {
      a = pi/2.;
      b = pi;
    }
    double fa = _Theta(Eps,L,Q,a,M,a);
    double c = (a+b)/2.;
    double fc = _Theta(Eps,L,Q,a,M,c);
    int N = 0;
    const int Nmax=1000;
    while( N < Nmax ) {
      c = (a+b)/2.; // New midpoint
      fc = _Theta(Eps,L,Q,a,M,c);
      if( c == 0 || (b-a)/2. < tolerance ) break;
      N++;
      
      if (sgn(fc) == sgn(fa)) {
        a = c;
        fa = fc;
      } else {
        b = c;
      }
    }
    // Check if solution found
    if( N < Nmax ) {
      thLimits[i] = c;
    } else {
      cerr << "WARNING: kerr_th_limits failed" << endl;
      thMin = std::numeric_limits<double>::quiet_NaN();
      thMax = std::numeric_limits<double>::quiet_NaN();
      return false;
    }
  }
  // Found the solution, input values!
  thMin = thLimits[0];
  thMax = thLimits[1];
  const double c = 0.5*(thMax+thMin);
  const double fc =_Theta(Eps,L,Q,a,M,c);
  if( fc < 0 ) {
    thMin = std::numeric_limits<double>::quiet_NaN();
    thMax = std::numeric_limits<double>::quiet_NaN();
    return false;
  }
  /////////
  return true;
}


bool kerr_r_limits( const double Eps,
                    const double L,
                    const double Q,
                    const double a,
                    const double M,
                    double & rMin,
                    double & rMax ) {
  const double a2=a*a;
  const double Eps2=Eps*Eps;
  const double L2  = L*L;
  //const double p=poly([0]);
  // Directly copied from Gyoto
  // Calculate roots of the effective potential:
  // c[0]u^4+... // See Mathematica notebook from 2.27.
  double results[MAXDEGREE] = {0};
  double results_im[MAXDEGREE] = {0};
  double coeffs[MDP1];
  int degree=4;
  // NOTE: These coefficients have been simplified with fullsimplify
  const double c4 = (Eps2-1.),                     // *r^4
               c3 = 2.*M,                          // *r^3
               c2 = -1.*L2-Q+a2*(Eps-1.)*(Eps+1.), // *r^2
               c1 = 2.*M*(Q+pow(L-a*Eps,2)),
               c0 = -1.*a2*Q;
  coeffs[0] = c4; // r^4
  coeffs[1] = c3; // r^3
  coeffs[2] = c2; // r^2
  if( coeffs[2] != coeffs[2] ) { cout << "Eps L Q " << Eps << " " << L << " " << Q << endl;; exit(1); }
  coeffs[3] = c1; //r
  coeffs[4] = c0; //c0
//  cout << "Coeffs: "; for( int i = 0; i < 5; ++i ) cout << coeffs[i] << " "; cout << endl;
  

  rpoly_ak1(coeffs, &degree, results, results_im);

  sort(results, results+4);
  // Save the values:
  bool resultsOk = true;
  for( int i = 0; i < 4; ++i ) { if( results_im[i] != 0 ) resultsOk = false; }
  {
    const double r = results[2]+FLOAT_ERR;//results[2]+0.00001*(results[3]-results[2]);
    const double r2= r*r, r3=r2*r, r4=r3*r;
    double soln= coeffs[0]*r4 + coeffs[1]*r3 + coeffs[2] * r2 + coeffs[3] * r + coeffs[4];
    if( soln < 0 ) resultsOk = false;
  }
  if( resultsOk ) {
    rMin = results[2];
    rMax = results[3];
  } else {
    rMin = std::numeric_limits<double>::quiet_NaN();
    rMax = std::numeric_limits<double>::quiet_NaN();
    return false;
  }
  const double rHorizon = M+sqrt(M*M-a*a);
  if(rMin <= rHorizon) return false; // The particle cant be inside horizon
  assert( results[2] <= results[3] );
  return true;
}


int I_th_integrand( unsigned ndim, const double *x, void *params_void,
		   unsigned fdim, double *fval ) {
  const double * params = (double*)params_void;
  const double th = *x;
  const double Eps = params[0];
  const double Q   = params[1];
  const double L   = params[2];
  const double a   = params[3];
  const double M   = params[4];
  const double sqrtTerm = _Theta(Eps,L,Q,a,M,th);
  if( sqrtTerm < 0 ) { *fval = 0; return 0;}
  *fval = sqrt(sqrtTerm);
  return 0;
}

double I_th(const double Eps,
          const double Q,
          const double L, 
          const double a, 
          const double M, 
          const double thmin, 
          const double thmax) {
  double params[5] = {Eps, Q, L, a, M};
  assert( thmax >= thmin );
  unsigned fdim = 1; unsigned dim = 1;
  size_t maxEval = 100000;
  double reqAbsError = 1.e-4;
  double reqRelError = 1.e-6;
  error_norm errorNorm;
  double result, error;
  integrand f = I_th_integrand;
  hcubature(fdim, f, params,
  dim, &thmin, &thmax,
  maxEval, reqAbsError, reqRelError, errorNorm,
  &result, &error
  );
  return result;
}

int I_r_integrand( unsigned ndim, const double *x, void *params_void,
		   unsigned fdim, double *fval ) {
  const double * params = (double*)params_void;
  const double Eps= params[0];
  const double Q   = params[1];
  const double L  = params[2];
  const double a  = params[3];
  const double M   = params[4];
  const double r = *x;
  const double r2  = r*r;
  const double a2  = a*a;
  const double Delta = r2 + a2 - 2.*G*M*r;
  const double R = _R(Eps,L,Q,a,M,r);
  if( R < 0 ) { *fval = 0; return 0;}
  *fval = sqrt(R) / Delta;
  return 0;
}


double I_r(const double Eps,
           const double L,
           const double Q,
           const double a,
           const double M,
           const double rmin,
           const double rmax) {
  assert( rmax >= rmin );
  double params[5] = {Eps, Q, L, a,M};
  unsigned fdim = 1; unsigned dim = 1;
  size_t maxEval = 100000;
  double reqAbsError = 1.e-6;
  double reqRelError = 1.e-6;
  error_norm errorNorm;
  double result, error;
  integrand f = I_r_integrand;
  hcubature(fdim, f, params,
  dim, &rmin, &rmax,
  maxEval, reqAbsError, reqRelError, errorNorm,
  &result, &error
  );
  return result;
}

double SysPrimeToTdot(const double pos[4], const double v[3],const double a,const double M) {
  double sum=0.,xpr[4];
  double g[4][4];
  int i,j;

  xpr[0]=1.; // dt/dt=1;
  for (i=0;i<3;++i) xpr[i+1]=v[i];
  gmunu(pos,a,M,g);
  for (i=0;i<4;++i) {
    for (j=0;j<4;++j) {
      sum+=g[i][j]*xpr[i]*xpr[j];
    }
  }
  if (sum>=0) {
    cout << "WARNING: Faster than speed of light at " << __FILE__ << " " << __LINE__ <<  endl;
    return 0.;
  }
  return pow(-sum, -0.5);
}

bool constants_to_pos( const long double Eps, const long double L, const long double Q, const long double a, const long double M, double x[4], double u[4], int dir ) {
  assert( dir == 1 || dir == -1 );
  double rmin, rmax, thmin, thmax;
  bool found_limits = kerr_r_limits(Eps,L,Q,a,M,rmin,rmax);
  if( found_limits == false ) return false;
  found_limits = kerr_th_limits(Eps,L,Q,a,M,thmin,thmax);
  if( found_limits == false ) return false;
  const double a2 = a*a;
  const double Eps2 = Eps*Eps;
  const double L2 = L*L;
  if( found_limits == false ) return false;
  // From claimedBest
  const double r=rmin+(rmax-rmin)/2.;
  const double r2= r*r;
  mt19937_64 gen(rand());
  double th, ph;
  uniform_real_distribution<> phi(0,2*pi);
  uniform_real_distribution<> costheta(-1.,1.);
  ph = phi(gen);
  while( true ) {
    th = acos(costheta(gen));
    if( th > thmin && th < thmax && _Theta(Eps,L,Q,a,M,th) > 0 ) break;
  }
  const double sinth2 = sin(th)*sin(th); const double sinth=sin(th);
  const double costh2 = cos(th)*cos(th);
  const double Delta=r2-2.*G*M*r+a2;
  const double Sigma=r2+a2*costh2;
  
  // Integrals of motion (8.5.10)
  double R    = _R(Eps,L,Q,a,M,r);
  if( R < 0 ) cout << "ERROR: R=" << R << " at " << __FILE__ << " " << __LINE__ << endl;
  double Omega= _Theta(Eps,L,Q,a,M,th);
  // Particle motion:
  {
    // Assign the direction randomly:
    double dtdlambda = (-a*(a*Eps*sinth2-L)+(r2+a2)/Delta * (Eps*(r2+a2)-L*a))/Sigma;
    dir = 2*(rand()%2)-1;
    double drdlambda = (double)dir*sqrt(R)/Sigma;
    dir = 2*(rand()%2)-1;
    double dthdlambda= (double)dir*sqrt(Omega)/Sigma;
    double dphdlambda= (-(a*Eps-L/sinth2)+a/Delta * (Eps*(r2+a2)-L*a))/Sigma;
    u[0] = dtdlambda;
    u[1] = drdlambda;
    u[2] = dthdlambda;
    u[3] = dphdlambda;
  }
  // Assign 4-position
  x[0] = 0;
  x[1] = r;
  x[2] = th;
  x[3] = 0;
  // Check for valid values
  for( int i = 0; i < 4; ++i ) {
    if( !(u[i] == u[i]) ) {cout << "u = " << u[0] << " " << u[1] << " " << u[2] << " " << u[3] << endl; cout << "Eps Lz Q: " << Eps << " " << L << " " << Q << endl; cout << "Omega: " << Omega << endl;}
    if( !(u[i] == u[i]) ) return false; // nan  check
  }
#ifndef NODEBUGGING
  // Check to make suire that these can be inversed
  double consts[4];
  constants(x, u, a, M, consts);
  if (abs(Eps-consts[0]) + abs(L  -consts[1]) + abs(Q  -consts[2]) > 1.e-6) return false;
  assert( abs(Eps-consts[0]) < 1.e-6 );
  assert( abs(L  -consts[1]) < 1.e-6 );
  assert( abs(Q  -consts[2]) < 1.e-6 );
#endif
  return true;
}

double nelder_mead_root( const gsl_vector * v, void * params_void ) {
  const double * params = (double*)params_void;
  const double Eps = gsl_vector_get(v,0);
  const double Q   = gsl_vector_get(v,1);
  //cout << "Nelder mead params: " << Eps << " " << Q << endl;
  // Get the known values of I_r and I_th
  const double L      = params[0];
  const double a_2    = params[1];
  const double M_2    = params[2];
  const double I_r_b  = params[3];
  const double I_th_b = params[4];
  
  //cout << "At nelder-mead: Trying out Eps, Q, L, a2, M2, I_r_b, I_th_b: "  << Eps << " " << Q << " " << L << " " << a_2 << " " << M_2 << " " << I_r_b << " " << I_th_b << endl;
  // Check to make sure energy doesnt go out of bounds:
  if( Eps > 1 || Q < 0 ) return numeric_limits<double>::max();

  // Get limits:
  double rmin, rmax, thmin, thmax;
  rmin=0; rmax=0;
  bool found_limits = kerr_r_limits(Eps,L, Q, a_2, M_2, rmin, rmax);
  if( found_limits == false || rmin < M_2 ){ return numeric_limits<double>::max();}
  found_limits = kerr_th_limits(Eps, L, Q, a_2, M_2, thmin, thmax);
  if( found_limits == false ) { return numeric_limits<double>::max(); }
  // If limits found, calculate value:
  const double I_r_a = I_r(Eps, L, Q, a_2, M_2, rmin, rmax);
  const double I_th_a= I_th(Eps, Q, L, a_2, M_2, thmin, thmax);
  //cout << "Managed to find the limits, Eps, Q, rmin,rmax, I_r_a, I_th_a: " << Eps << " " << Q  << " " << rmin << " " << rmax << " " << I_r_a << " " << I_th_a << endl;
  //cout << "Value: " <<  abs(I_r_b - I_r_a) + abs(I_th_b - I_th_a) << endl;
  return abs(I_r_b - I_r_a)/abs(I_r_b) + abs(I_th_b - I_th_a)/abs(I_th_b);
}

bool startingValues(const double Epsb, 
        const double L, 
        const double Qb, 
        const double a_1,
        const double M_1,
        const double a_2, 
        const double M_2, 
        double & EpsStart, 
        double & QStart) {
  double rmin, rmax, thmin, thmax;
  assert( kerr_r_limits(Epsb,L,Qb,a_1,M_1,rmin,rmax) );
  random_device rd;
  mt19937 mt(rd());
  uniform_real_distribution<double> distE (0.95,1.0);
  uniform_real_distribution<double> distQ(0,(60.*M_2)*(60.*M_2));

  // Find a starting value by brute force
  const int imax = 10000000;
  int i;
  for( i = 0; i < imax; ++i ) {
    EpsStart = distE(mt);
    QStart   = distQ(mt);
    if( kerr_r_limits(EpsStart, L, QStart, a_2, M_2, rmin, rmax) && kerr_th_limits(EpsStart, L, QStart, a_2, M_2, thmin, thmax ) && rmin < RHARDCUTOFF*M_2 ) break;
  }
  assert( i < imax );

  if( i >= imax ) { 
#ifndef NDEBUG
    cerr << "WARNING: startingValues failed to find starting value!" <<endl; 
#endif
    return false;
    
  } // Failed to find starting value
  return true;
}





//using namespace libcmaes;
//
//namespace cma_solver {
//
//int solve( const double EpsStart, const double QStart, const double params[5], double & Epsa, double & Qa, double & result ) {
//  // Declare function ( I know this looks ugly )
//  FitFunc cma_root = [params, QStart](const double *x, const int N)
//  {
//    const double Eps = x[0]/QStart;
//    const double Q   = x[1];
//    if( Eps >= 1 || Q < 0 ) return numeric_limits<double>::max();
//    // Get the known values of I_r and I_th
//    const double L      = params[0];
//    const double a_2    = params[1];
//    const double M_2    = params[2];
//    const double I_r_b  = params[3];
//    const double I_th_b = params[4];
//
//    // Check to make sure energy doesnt go out of bounds:
//    if( Eps < 0 || Eps > 1 || Q < 0 ) return numeric_limits<double>::max();
//
//    // Get limits:
//    double rmin, rmax, thmin, thmax;
//    bool found_limits = kerr_r_limits(Eps,L, Q, a_2, M_2, rmin, rmax);
//    if( found_limits == false ){ return numeric_limits<double>::max();}
//    found_limits = kerr_th_limits(Eps, L, Q, a_2, M_2, thmin, thmax);
//    if( found_limits == false ) { return numeric_limits<double>::max(); }
//    // If limits found, calculate value:
//    const double I_r_a = I_r(Eps, L, Q, a_2, M_2, rmin, rmax);
//    const double I_th_a= I_th(Eps, Q, L, a_2, M_2, thmin, thmax);
//    //cout << "Managed to find the limits, Eps, Q, rmin,rmax, I_r_a, I_th_a: " << Eps << " " << Q  << " " << rmin << " " << rmax << " " << I_r_a << " " << I_th_a << endl;
//    //cout << "Value: " <<  abs(I_r_b - I_r_a) + abs(I_th_b - I_th_a) << endl;
//    return abs(I_r_b - I_r_a)/abs(I_r_b) + abs(I_th_b - I_th_a)/abs(I_th_b);
//  };
//  const int dim = 2; // problem dimensions.
//  std::vector<double> x0; x0.resize(dim); x0[0] = QStart*EpsStart; x0[1] = QStart;
//  double sigma = 0.01*QStart;
//  //int lambda = 100; // offsprings at each generation.
//  CMAParameters<> cmaparams(x0,sigma);
//  //cmaparams.set_algo(BIPOP_CMAES);
//  CMASolutions cmasols = cmaes<>(cma_root,cmaparams);
//  Candidate soln = cmasols.best_candidate();
//  x0 = soln.get_x();
//  Epsa = x0[0]/QStart;
//  Qa = x0[1];
//  result = soln.get_fvalue();
//  if( abs(result) < SUCCESSFULTOLERANCE ) return SUCCESSFUL;
//  return UNSUCCESSFUL;
//}
//}

namespace NonLinear {
  
double previous_y0[60];
double previous_y1[60];



int
print_state (size_t iter, gsl_multiroot_fsolver * s)
{
  printf ("iter = %3u x = % .3f % .3f "
          "f(x) = % .3e % .3e\n",
          iter,
          gsl_vector_get (s->x, 0), 
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->f, 0), 
          gsl_vector_get (s->f, 1));
}

int
nonlinear_root(const gsl_vector * x, void *params_void, 
              gsl_vector * f)
{
  const int tid = omp_get_thread_num();
  const double * params = (double*)params_void;

  // Get the known values of I_r and I_th
  const double L      = params[0];
  const double a_2    = params[1];
  const double M_2    = params[2];
  const double I_r_b  = params[3];
  const double I_th_b = params[4];

  const double Eps = gsl_vector_get (x, 0);
  const double Q   = gsl_vector_get (x, 1);
  if( Eps != Eps || Q != Q ) { return GSL_FAILURE; } // Make sure there are no NaNs
  //cout << "Eps Q: " << Eps << " " << Q << endl;
  // Check this out for non-linear optimizations: http://www.aip.de/groups/soe/local/numres/bookcpdf/c9-7.pdf
  if( Eps >= 1 || Q < 0 ) { 
    gsl_vector_set (f, 0, previous_y0[tid]); 
    gsl_vector_set (f, 1, previous_y1[tid]); 
    return GSL_SUCCESS;
  }
  // Check to make sure energy doesnt go out of bounds:
  if( Eps < 0 || Eps > 1 || Q < 0 ) { 
    gsl_vector_set (f, 0, previous_y0[tid]); 
    gsl_vector_set (f, 1, previous_y1[tid]); 
    return GSL_SUCCESS;
  }

  // Get limits:
  double rmin, rmax, thmin, thmax;
  bool found_limits = kerr_r_limits(Eps,L, Q, a_2, M_2, rmin, rmax);
  if( found_limits == false ){ 
    gsl_vector_set (f, 0, previous_y0[tid]); 
    gsl_vector_set (f, 1, previous_y1[tid]); 
    return GSL_SUCCESS;
  }
  found_limits = kerr_th_limits(Eps, L, Q, a_2, M_2, thmin, thmax);
  if( found_limits == false ) { 
    gsl_vector_set (f, 0, previous_y0[tid]); 
    gsl_vector_set (f, 1, previous_y1[tid]); 
    return GSL_SUCCESS;
  }
  // If limits found, calculate value:
  const double I_r_a = I_r(Eps, L, Q, a_2, M_2, rmin, rmax);
  const double I_th_a= I_th(Eps, Q, L, a_2, M_2, thmin, thmax);

  const double y0 = I_r_a - I_r_b;
  const double y1 = I_th_a - I_th_b;

  gsl_vector_set (f, 0, y0);
  gsl_vector_set (f, 1, y1);
  previous_y0[tid] = y0;
  previous_y1[tid] = y1;

  return GSL_SUCCESS;
}


int solve( const double EpsStart, const double QStart, double params[5], double & Epsa, double & Qa, double & result ) {
  double rmin, rmax, thmin, thmax;
  // Get the known values of I_r and I_th
  const double L      = params[0];
  const double a_2    = params[1];
  const double M_2    = params[2];
  const double I_r_b  = params[3];
  const double I_th_b = params[4];
 
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;

  int status;
  size_t i, iter = 0;

  const size_t n = 2;
  
  gsl_multiroot_function f = {&nonlinear_root, n, &params[0]};

  double x_init[2] = {EpsStart, QStart};
  gsl_vector *x = gsl_vector_alloc (n);

  gsl_vector_set (x, 0, x_init[0]);
  gsl_vector_set (x, 1, x_init[1]);

  T = gsl_multiroot_fsolver_hybrids;
  s = gsl_multiroot_fsolver_alloc (T, 2);
  gsl_multiroot_fsolver_set (s, &f, x);

  print_state (iter, s);

  do
    {
      iter++;
      status = gsl_multiroot_fsolver_iterate (s);

      //print_state (iter, s);

      if (status)   /* check if solver is stuck */
        break;

      status = 
        gsl_multiroot_test_residual (s->f, 1e-7);
    }
  while (status == GSL_CONTINUE && iter < 1000);

  const double y0 = gsl_vector_get(s->f,0);
  const double y1 = gsl_vector_get(s->f,1);
  result = abs(y0)/I_r_b+abs(y1)/I_th_b;
  Epsa = gsl_vector_get(s->x,0);
  Qa = gsl_vector_get(s->x,1);

  printf ("status = %s %.16e %.16e %.16e \n", gsl_strerror (status), result, Epsa, Qa);

  gsl_multiroot_fsolver_free (s);
  gsl_vector_free (x);
  
  if( result < SUCCESSFULTOLERANCE ) return SUCCESSFUL;
  return UNSUCCESSFUL;
}
}

namespace HomeBrewNewtonMethod { 
using namespace Eigen;
namespace DD {
const double dIrdE(const double Eps, const double L, const double Q, const double a, const double M) {
  const double h = sqrt(DBL_EPSILON);
  double div = 2.*h;
  double Eps2 = Eps + h;
  double Q2   = Q   + h;
  double Eps1 = Eps - h;
  double Q1   = Q   - h;
  if( Eps2 >= 1 ) {Eps2 = Eps; Q2 = Q; div = h;}
  double rmin1, rmin2, rmax1, rmax2;
  bool found_limits1 = kerr_r_limits(Eps1,L, Q, a, M, rmin1, rmax1);
  bool found_limits2 = kerr_r_limits(Eps2,L, Q, a, M, rmin2, rmax2);
  // Make sure they're not out of bounds!
  if( found_limits1 == false && found_limits2 == false ) { return std::numeric_limits<double>::quiet_NaN(); }
  if( found_limits1 == false ) {Eps1 = Eps; div = h;}
  if( found_limits2 == false ) {Eps2 = Eps; div = h;}

  // Calculate the derivative
  const double Ir1 = I_r(Eps1, L, Q, a, M, rmin1, rmax1 );
  const double Ir2 = I_r(Eps2, L, Q, a, M, rmin2, rmax2 );
  return (Ir2-Ir1)/div;
}

const double dIrdQ(const double Eps, const double L, const double Q, const double a, const double M) {
  const double h = sqrt(DBL_EPSILON);
  double div = 2.*h;
  double Eps2 = Eps + h;
  double Q2   = Q   + h;
  double Eps1 = Eps - h;
  double Q1   = Q   - h;
  if( Q1 < 0 ) {Q1 = Q; div = h;}
  double rmin1, rmin2, rmax1, rmax2;
  bool found_limits1 = kerr_r_limits(Eps,L, Q1, a, M, rmin1, rmax1);
  bool found_limits2 = kerr_r_limits(Eps,L, Q2, a, M, rmin2, rmax2);
  // Make sure they're not out of bounds!
  if( found_limits1 == false && found_limits2 == false ) { return std::numeric_limits<double>::quiet_NaN(); }
  if( found_limits1 == false ) {Q1 = Q; div = h;}
  if( found_limits2 == false ) {Q2 = Q; div = h;}
  // Calculate the derivative
  const double Ir1 = I_r(Eps, L, Q1, a, M, rmin1, rmax1 );
  const double Ir2 = I_r(Eps, L, Q1, a, M, rmin2, rmax2 );
  return (Ir2-Ir1)/div;
}
const double dIthdQ(const double Eps, const double L, const double Q, const double a, const double M) {
  const double h = sqrt(DBL_EPSILON);
  double div = 2.*h;
  double Eps2 = Eps + h;
  double Q2   = Q   + h;
  double Eps1 = Eps - h;
  double Q1   = Q   - h;
  if( Q1 < 0 ) {Q1 = Q; div = h;}
  double thmin1, thmin2, thmax1, thmax2;
  bool found_limits1 = kerr_th_limits(Eps,L, Q1, a, M, thmin1, thmax1);
  bool found_limits2 = kerr_th_limits(Eps,L, Q2, a, M, thmin2, thmax2);
  // Make sure they're not out of bounds!
  if( found_limits1 == false && found_limits2 == false ) { return std::numeric_limits<double>::quiet_NaN(); }
  if( found_limits1 == false ) {Q1 = Q; div = h;}
  if( found_limits2 == false ) {Q2 = Q; div = h;}
  // Calculate the derivative
  const double Ir1 = I_th(Eps, Q1, L, a, M, thmin1, thmax1 );
  const double Ir2 = I_r(Eps, Q1, L, a, M, thmin2, thmax2 );
  return (Ir2-Ir1)/div;
}
const double dIthdE(const double Eps, const double L, const double Q, const double a, const double M) {
  const double h = sqrt(DBL_EPSILON);
  double div = 2.*h;
  double Eps2 = Eps + h;
  double Q2   = Q   + h;
  double Eps1 = Eps - h;
  double Q1   = Q   - h;
  if( Eps2 >= 1 ) {Eps2 = Eps; Q2 = Q; div = h;}
  double rmin1, rmin2, rmax1, rmax2;
  bool found_limits1 = kerr_th_limits(Eps1,L, Q, a, M, rmin1, rmax1);
  bool found_limits2 = kerr_th_limits(Eps2,L, Q, a, M, rmin2, rmax2);
  // Make sure they're not out of bounds!
  if( found_limits1 == false && found_limits2 == false ) { return std::numeric_limits<double>::quiet_NaN(); }
  if( found_limits1 == false ) {Eps1 = Eps; div = h;}
  if( found_limits2 == false ) {Eps2 = Eps; div = h;}

  // Calculate the derivative
  const double Ir1 = I_th(Eps1, Q, L, a, M, rmin1, rmax1 );
  const double Ir2 = I_th(Eps2, Q, L, a, M, rmin2, rmax2 );
  return (Ir2-Ir1)/div;
}
}

int solve( const double EpsStart, const double QStart, double params[5], double & Epsa, double & Qa, double & result ) {
  // Get the known values of I_r and I_th
  const double L      = params[0];
  const double a_2    = params[1];
  const double M_2    = params[2];
  const double I_r_b  = params[3];
  const double I_th_b = params[4];
  
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;

  int status;
  size_t i, iter = 0;

  const size_t dim = 2;
  
  // Set the initial Eps and Q values
  double Eps = EpsStart;
  double Q = QStart;
  double f;
  const double max_iterations = 1000;
  for( int i = 0; i < max_iterations; ++i ) {
    // Get limits:
    double rmin, rmax, thmin, thmax;
    bool found_limits = kerr_r_limits(Eps,L, Q, a_2, M_2, rmin, rmax);
    if( found_limits == false ){  }
    found_limits = kerr_th_limits(Eps, L, Q, a_2, M_2, thmin, thmax);
    if( found_limits == false ) {  }
    // If limits found, calculate value:
    const double I_r_a = I_r(Eps, L, Q, a_2, M_2, rmin, rmax);
    const double I_th_a= I_th(Eps, Q, L, a_2, M_2, thmin, thmax);
    const double y0 = I_r_b - I_r_a;
    const double y1 = I_th_b- I_th_a;
    Eigen::Matrix<double,2,1> F(y0,y1);
    Eigen::Matrix<double,2,2> J;
    {
      const double dF00 = DD::dIrdE(Eps,L,Q,a_2,M_2);
      const double dF01 = DD::dIrdQ(Eps,L,Q,a_2,M_2);
      const double dF10 = DD::dIthdE(Eps,L,Q,a_2,M_2);
      const double dF11 = DD::dIthdQ(Eps,L,Q,a_2,M_2);
      if( dF00 != dF00 || dF01 != dF01 || dF10 != dF10 || dF11 != dF11 ) return UNSUCCESSFUL; // Make sure not nan!
      J << dF00, dF01,
           dF10, dF11;
    }
    Eigen::HouseholderQR<Matrix2d> qr(J);
    Eigen::Matrix<double,2,1> dx = -1.*qr.solve(F);
    // backtrack: http://www.aip.de/groups/soe/local/numres/bookcpdf/c9-7.pdfa
    double lambda = 1; const double descent = 0.9;
    const double fold = 0.5*(pow(I_r_b - I_r_a,2)+pow(I_th_b-I_th_a,2));
    double Epsnew, Qnew;
    while( true ) {
      // Use Newton's method to find the criteria
      Epsnew = Eps + lambda*dx(0);
      Qnew = Q + lambda*dx(1);
      // Get limits:
      double rmin, rmax, thmin, thmax;
      bool found_limits = kerr_r_limits(Epsnew,L, Qnew, a_2, M_2, rmin, rmax);
      if( found_limits == false ){ lambda *= descent; continue; }
      found_limits = kerr_th_limits(Epsnew, L, Qnew, a_2, M_2, thmin, thmax);
      if( found_limits == false ) { lambda *= descent; continue; }
      if( Epsnew >= 1 || Q < 0 ) { lambda *= descent; continue; }
      // If limits found, calculate value:
      const double I_r_a_new = I_r(Epsnew, L, Qnew, a_2, M_2, rmin, rmax);
      const double I_th_a_new= I_th(Epsnew, Qnew, L, a_2, M_2, thmin, thmax);
      f = 0.5*(pow(I_r_b - I_r_a_new,2)+pow(I_th_b-I_th_a_new,2)); // TODO: Calculate these!
      //cout << f << " " << fold << endl;
      if( f <= fold ) break; // Ok value
      lambda *= descent;
    }
    // Update values:
    Eps = Epsnew;
    Q   = Qnew;
    if( f < 1.0e-7 ) break; // Check value
  }
  
  result = f;
  
  if( result < 1.e-6 ) return SUCCESSFUL;
  return UNSUCCESSFUL;
}
}

// namespace Annealing {
// 
//   
// double root_solver( void * xp ) {
//   double * x = ((double*)xp);
//   
//   return 0;
// }
// int solve( const double EpsStart, const double QStart, double params[5], double & Epsa, double & Qa, double & result ) {
//   const gsl_rng_type * T;
//   gsl_rng * r;
// 
//   double x_initial = 15.5;
// 
//   gsl_rng_env_setup();
// 
//   T = gsl_rng_default;
//   r = gsl_rng_alloc(T);
// 
//   gsl_siman_solve(r, &x_initial, root_solver, S1, M1, P1,
//                   NULL, NULL, NULL, 
//                   sizeof(double), params);
// 
//   gsl_rng_free (r);
//   return 0;
// }
// 
// }



namespace SimulatedAnnealing {
double drand() {return (double)rand()/(double)RAND_MAX;}

double P(double E1, double E2, double T) {
  if( E2<E1 ) return 1;
  else        return exp(-1.*(E2-E1)/T);
}

int solve( const double EpsStart, const double QStart, double params[5], double & Epsa, double & Qa, double & result ) {
  double rmin, rmax;
  // Get the known values of I_r and I_th
  const double L      = params[0];
  const double a_2    = params[1];
  const double M_2    = params[2];
  const double I_r_b  = params[3];
  const double I_th_b = params[4];
  gsl_vector *v, *vnew;
  v =  gsl_vector_alloc(2);
  vnew=gsl_vector_alloc(2);

  gsl_vector_set(v, EpsStart, 0);
  gsl_vector_set(v, QStart,   1);
  result = nelder_mead_root( v, (void*)params );

  const int kmax = 1000000;
  for( int k = 0; k < kmax; ++k ) {
    double Eps = gsl_vector_get(v, 0);
    double Q   = gsl_vector_get(v, 1);
    double EpsNew= Eps+(2.*drand()-1.)*0.00001;
    double QNew  = Q  +(2.*drand()-1.)*1;
    bool found_limits = kerr_r_limits(Eps, L, Q,a_2,M_2,rmin,rmax);
    if( found_limits == false ) { continue;}
    double T = (1.-(double)k/(double)kmax)*10;
    // Set new vector
    gsl_vector_set(vnew, EpsNew, 0);
    gsl_vector_set(vnew, QNew,   1);
    // Check if we accept state:
    double Enew = nelder_mead_root( vnew, (void*)params );
    double E    = nelder_mead_root( v, (void*)params );
    if( P(E, Enew, T) >= drand() ) {
      gsl_vector_set(v, EpsNew, 0);
      gsl_vector_set(v, QNew,   0);
    }
  }

  // Show result
  result = nelder_mead_root( v, (void*)params );

  if ( result < SUCCESSFULTOLERANCE ) {
    return SUCCESSFUL;
  }
  return UNSUCCESSFUL;
}

}

namespace NelderMead {
int solve( const double EpsStart, const double QStart, double params[5], double & Epsa, double & Qa, double & result ) {
  gsl_set_error_handler_off();
  // Get the known values of I_r and I_th
  const double L      = params[0];
  const double a_2    = params[1];
  const double M_2    = params[2];
  const double I_r_b  = params[3];
  const double I_th_b = params[4];
  // Initialize
  const gsl_multimin_fminimizer_type *T = 
    gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *X;
  gsl_multimin_function minex_func;

  size_t iter = 0;
  int status;
  double size;
  
  X = gsl_vector_alloc (2);
  gsl_vector_set (X, 0, EpsStart);
  gsl_vector_set (X, 1, QStart);

  /* Set initial step sizes to 1 */
  ss = gsl_vector_alloc (2);
  gsl_vector_set_all (ss, 1.0);

  /* Initialize method and iterate */
  minex_func.n = 2;
  minex_func.f = nelder_mead_root;
  minex_func.params = params;

  s = gsl_multimin_fminimizer_alloc (T, 2);
  gsl_multimin_fminimizer_set (s, &minex_func, X, ss);

  do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);
      
      if (status) 
        break;

      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, 1e-6);

      if (status == GSL_SUCCESS)
        {
        }
    }
  while (status == GSL_CONTINUE && iter < 10000);
  
  result = s->fval;
  Epsa = gsl_vector_get(s->x,0);
  Qa   = gsl_vector_get(s->x,1);
  
  gsl_vector_free(X);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);
  if( status != GSL_SUCCESS){ cout << "GSL Integration failed with Epsa Qa: " << Epsa << " " << Qa << endl; return UNSUCCESSFUL;}
  if ( result < SUCCESSFULTOLERANCE ) {
    return SUCCESSFUL;
  }
  return UNSUCCESSFUL;
}
}

int grow_black_hole(const double consts_1[4],
                    const double a_1,
                    const double a_2,
                    const double M_1,
                    const double M_2,
                    const Solver solver,
                    double consts_2[4],
                    double & result) {
#ifndef DONTCHECKCONDITIONS
  assert( check_conditions( M_1, M_2, a_1, a_2 ) );
#endif
  // Get constants:
  const double L  = consts_1[1];
  const double Epsb=consts_1[0];
  const double Qb = consts_1[2];
  if( Epsb >= 1 || Qb < 0 ) {cout <<" Bad eps: " << Epsb << endl; return UNSUCCESSFUL;}
  // Get limits:
  double rmin, rmax, thmin, thmax;
  bool found_limits = kerr_r_limits(Epsb, L, Qb,a_1,M_1,rmin,rmax);
  if( found_limits == false ) { cout << "Bad r limits: " << rmin << " " << rmax << endl << endl; return UNSUCCESSFUL;}
  found_limits = kerr_th_limits(Epsb, L, Qb, a_1,M_1,thmin,thmax);
  if( found_limits == false ) { cout << "Bad th limits: " << thmin << " " << thmax << endl << endl; return UNSUCCESSFUL; }
  // Adiabatic invariants:
  const double I_r_b = I_r(Epsb,L,Qb,a_1,M_1, rmin, rmax);
  const double I_th_b= I_th(Epsb, Qb,L, a_1,M_1, thmin, thmax);

  /* Starting point */
  double EpsStart, QStart;
  // TODO: Make the initial guess non-adiabatic growth
  //assert( startingValues(Epsb, L, Qb, a_1, M_1, a_2, M_2, EpsStart, QStart) );
  if( startingValues(Epsb, L, Qb, a_1, M_1, a_2, M_2, EpsStart, QStart) == false ) {cout << "Bad starting value" << endl; return UNSUCCESSFUL;}
  {
    double rminb,rmaxb, thminb, thmaxb;
    kerr_r_limits(EpsStart,L,QStart,a_2,M_2,rminb,rmaxb);
    kerr_th_limits(EpsStart,L,QStart,a_2,M_2,thminb,thmaxb);
    assert( kerr_r_limits(EpsStart,L, QStart, a_2, M_2, rminb, rmaxb) && kerr_th_limits(EpsStart, L, QStart, a_2, M_2, thminb, thmaxb ) );
  }


  /* Initialize method and iterate */
  double params[5] = {L, a_2, M_2, I_r_b, I_th_b};
  
  // Solve the equation:
  double Epsa, Qa;
  int status;
  switch(solver) {
/*
    case CMA : {
      for( int i = 0; i < 40; ++i ) {
        status = cma_solver::solve( EpsStart, QStart, params, Epsa, Qa, result ); 
        if( status == SUCCESSFUL ) break;
        EpsStart = Epsa;
        QStart   = Qa;
      }
      if(status == UNSUCCESSFUL) {
        EpsStart = Epsa; QStart = Qa;
        status = NonLinear::solve(EpsStart, QStart, params, Epsa, Qa, result );
      }
    } break;
*/
    case nelderMead : status = NelderMead::solve(EpsStart, QStart, params, Epsa, Qa, result); break;
    case homebrew : status = HomeBrewNewtonMethod::solve(EpsStart, QStart, params, Epsa, Qa, result); break;
    case hybrid : {
      status = NelderMead::solve(EpsStart, QStart, params, Epsa, Qa, result);
      if(status == UNSUCCESSFUL) {
        EpsStart = Epsa; QStart = Qa;
        status = NonLinear::solve(EpsStart, QStart, params, Epsa, Qa, result );
        // TODO:
        // If the solution is *still* not found, go for random shooting:
        if( status == UNSUCCESSFUL ) {
          
        }
      }
    } break;
    case nonlinear : {
      cout << "Not implemented" << endl;
      status = NonLinear::solve(EpsStart, QStart, params, Epsa, Qa, result );
      status = UNSUCCESSFUL;
    }
    break;
  }
  // If the result is still unsuccessful, just return UNSUCCESSFUL
  if( status == UNSUCCESSFUL ) return UNSUCCESSFUL;

  // Get the new limits:
  found_limits = kerr_r_limits(Epsa,L, Qa, a_2, M_2, rmin, rmax);
  assert( found_limits == true);

  // Return the limits
  consts_2[0] = Epsa;
  consts_2[1] = L;
  consts_2[2] = Qa;
  const double Ka = Qa + pow(Epsa*a_2-L,2);
  consts_2[3] = Ka;
  const double I_r_a = I_r(Epsa,L,Qa,a_2,M_2, rmin, rmax);
  const double I_th_a= I_th(Epsa, Qa,L, a_2,M_2, thmin, thmax);
  if ( result < SUCCESSFULTOLERANCE ) {
    if( abs(I_r_a - I_r_b)/abs(I_th_b) + abs(I_th_a - I_th_b)/abs(I_th_b) > SUCCESSFULTOLERANCE ) cout << "Note: Found variables but the result is: " << abs(I_r_a - I_r_b)/abs(I_th_b) + abs(I_th_a - I_th_b)/abs(I_th_b) << " " << result << endl;
    return SUCCESSFUL;
  }
  return UNSUCCESSFUL;
}

} // End namespace helper



