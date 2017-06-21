#ifndef _GROWBLACKHOLE_H
#define _GROWBLACKHOLE_H
#include "sampling.h"
#include "lookupTable.h"
#include <random>

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
                   std::mt19937_64 & gen,
                   std::normal_distribution <double > & normal_dist,
                   LookupTable & lookup,
                   int & Nsampled,
                   int & Nin,
                   std::vector<double> & EpsVectorNew,
                   std::vector<double> & LzVectorNew,
                   std::vector<double> & QVectorNew,
                   std::vector<double> & EpsVector,
                   std::vector<double> & LzVector,
                   std::vector<double> & QVector,
                   std::vector<double> & MStartVector,
                   std::vector<double> & aStartVector,
                   const bool sampleUnsuccessfulGrowths=false
);
#endif