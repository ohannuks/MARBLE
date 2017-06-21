#ifndef _LOOKUPTABLE_H
#define _LOOKUPTABLE_H
#include <vector>
#include <utility>
#include <unordered_map>

class LookupTable {
public:
  std::vector<double> EpsVector;
  std::vector<double> LzVector;
  std::vector<double> QVector;
  std::vector<std::pair<double,int>> Sr ;
  std::vector<std::pair<double,int>> Sth;
  std::vector<std::pair<double,int>> Lz;
  int nParticles;
  double rMinLimit;
  double a, M;
  int GrowBlackHole(const double Sr_1, const double Sth_1, const double Lz_1, const double M_1, const double a_1, const double M_2, const double a_2, double consts[4]);
  double maximumSr();
  double minimumSr();
  double maximumSth();
  double minimumSth();
  double maximumLz();
  double minimumLz();

};

int findKey(const double x, 
            const std::vector<std::pair<double,int>> & values,
            std::unordered_map<int,int> & keys );
bool read_lookup( LookupTable & lookup,
                  const char * filename
                 );

#endif