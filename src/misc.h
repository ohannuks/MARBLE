#include <string>
#include <vector>
#include <fstream>
#include <inttypes.h>
void save_results( 
                   const std::vector<double> & EpsVector,
                   const std::vector<double> & LzVector,
                   const std::vector<double> & QVector, 
                   const std::vector<double> & EpsVectorNew,
                   const std::vector<double> & LzVectorNew,
                   const std::vector<double> & QVectorNew, 
                   const std::vector<double> & MStartVector,
                   const std::vector<double> & aStartVector,
                   const long long Ntotal, 
                   const long long Nin,
                   const double a_1,
                   const double M_1,
                   const double a_2,
                   const double M_2,
                   const std::string & filename);

