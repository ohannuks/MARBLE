#include "misc.h"
#include <cassert>
#include <inttypes.h>
#define __STDC_FORMAT_MACROS
using namespace std;

void save_results( 
                   const vector<double> & EpsVector,
                   const vector<double> & LzVector,
                   const vector<double> & QVector, 
                   const vector<double> & EpsVectorNew,
                   const vector<double> & LzVectorNew,
                   const vector<double> & QVectorNew, 
                   const vector<double> & MStartVector,
                   const vector<double> & aStartVector,
                   const long long Ntotal, 
                   const long long Nin,
                   const double a_1,
                   const double M_1,
                   const double a_2,
                   const double M_2,
                   const string & filename) {
  assert( EpsVector.size() == QVector.size());
  FILE * f = fopen(filename.c_str(), "w");
  // Add attributes
  fprintf(f,"%lld \n", Ntotal);
  fprintf(f,"%lld \n", Nin);
  fprintf(f,"%.16e \n", a_1);
  fprintf(f,"%.16e \n", M_1);
  fprintf(f,"%.16e \n", a_2);
  fprintf(f,"%.16e \n", M_2);
  for( int i = 0; i < EpsVector.size(); ++i ) {
    fprintf(f,"%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",EpsVector[i], LzVector[i], QVector[i], EpsVectorNew[i], LzVectorNew[i], QVectorNew[i], MStartVector[i], aStartVector[i]);
  }
  fclose(f);
}



