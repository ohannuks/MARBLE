#ifndef _HELPER_FUNCS_
#define _HELPER_FUNCS_
#include <gsl/gsl_multimin.h>
#include <iostream>
#include "definition.h"


enum Solver{ CMA, nelderMead, hybrid, nonlinear, homebrew, bruteforce };



namespace helper {

/*
 *  Find u^0 from u^\mu u_\mu = -1 for converting between newtonian and non-newtonian approach
 */
void findVt( const double x[4], const double v[3], const double a, const double M, double u[4] );

void gmunu( const double x[4],
            const double a,
            const double M,
            double g[4][4] );
  /**
   * Gives the tetrad (vierbein) according to https://arxiv.org/pdf/1105.0094.pdf
   * 
   * \param pos Particle's 4-position
   * \param tetrad Returns the tetrad
   **/
  void ZAMOTetrad( double const x[4], const double a, const double M,
		   double tetrad[4][4], bool inverse=false );

/**
 * Returns Eps, L, Q, K
 **/
void constants(const double pos[4],
               const double u[4],
               const double a,
               const double M,
               double consts[4]);

bool constants_to_pos( const long double Eps, const long double L, const long double Q, const long double a, const long double M, double x[4], double u[4], int dir = 1 );

bool kerr_th_limits(const double Eps,
                    const double L,
                    const double Q,
                    const double a,
                    const double M,
                    double & thMin,
                    double & thMax );

double SysPrimeToTdot(const double pos[4], const double v[3],const double a,const double M);

bool kerr_r_limits( const double Eps,
                    const double L,
                    const double Q,
                    const double a,
                    const double M,
                    double & rMin,
                    double & rMax );

double I_th_integrand( double th, void * params_void );

double I_th(const double Eps,
          const double Q,
          const double L, 
          const double a, 
          const double M, 
          const double thmin, 
          const double thmax);

double I_r_integrand( const double r, void * params_void );

double I_r(const double Eps,
           const double L,
           const double Q,
           const double a,
           const double M,
           const double rmin,
           const double rmax);


double nelder_mead_root( const gsl_vector * v, void * params_void );

int  nelder_mead( const double x[4],
                  const double u[4],
                  const double a_1,
                  const double a_2,
                  const double M_1,
                  const double M_2,
                  double xnew[4],
                  double unew[4],
                  double & result);

int grow_black_hole(const double consts_1[4],
                    const double a_1,
                    const double a_2,
                    const double M_1,
                    const double M_2,
                    const Solver solver,
                    double consts_2[4],
                    double & result);
}












#endif
