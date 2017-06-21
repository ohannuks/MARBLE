#ifndef _SAMPLING_H_
#define _SAMPLING_H_
// Gyoto
#include "cstdlib"
#include <iostream>
#include <chrono> 
#include <ctime>
#include <cstdio>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>
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
#include "cubature/cubature.h"
#include <random>
#include "rpoly_ak1.h"
#include <algorithm>
#include <gsl/gsl_integration.h>
#include "definition.h"

namespace sampling {
  
enum Profile {constant, r1, r2, r3}; // Constant density profile, r^-1 density profile and r^-2, et cetera
  
/**
 * Maxwellian velocity sampling (note: we assume a locally isotropic distribution function )
 * 
 * \param normal_dist Normal gaussian distribution with mean value = 0 and stddev equal to velocity dispersion
 * \param gen A random number generator
 * \param mean_v Our mean velocity so the output will be mean_v+variance
 * \param r The position from which to shoot the particle (usually retrieved via r_sampling )
 * \param v The velocity is output here
 * \returns SUCCESSFUL when successful
 **/
int maxwellian_sampling( std::normal_distribution<double> & normal_dist, std::mt19937_64 & gen, const double mean_v[3], const double x[4], double v[3] ); 
  
/**
 * Samples the r paramter as a function of the INFLUENCE radius
 * \param rMinInfluence The influence radius when the black hole is an infant
 * \param rMaxInfluence The influence radius when the black hole is grown up
 * \param gen Random number generator (mersenne twister)
 * \param r The radius r
 * \returns SUCCESSFUL when successful
 **/
int r_sampling( const double rMinInfluence, const double rMaxInfluence, std::mt19937_64 & gen, double x[4], Profile profile, const double rMaxSample );
  




}



#endif












