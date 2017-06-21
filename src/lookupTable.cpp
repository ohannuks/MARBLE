#include "lookupTable.h"
#include <algorithm>
#include <limits>
#include "helperFunctions.h"
#include <cassert>
#include "definition.h"

using namespace std;

int findKey(const double x, 
            const vector<pair<double,int>> & values,
            unordered_map<int,int> & keys ) {
  vector<pair<double,int>>::const_iterator low,up;
  const double errTol = 0.02; // Percentage error
  low=std::lower_bound (values.begin(), values.end(), make_pair<double,int>(x-errTol*x,1) );
  up =std::upper_bound (values.begin(), values.end(), make_pair<double,int>(x+errTol*x,1) );
  for( int i = low-values.begin(); i < up-values.begin(); ++i ) {
    const int key = values[i].second;
    if( keys.find(key) == keys.end() ) keys[key] =  1; //Couldnt find key
    else                               keys[key] += 1;
    if( keys[key] == 3 ) return key; // Found three matches within error bounds
  }
  cout << "Lower and upper bound " << low-values.begin() << " " << up-values.begin() << " " << values.size() << endl;
  return -1;
}

double LookupTable::maximumSr()
{
  return Sr[Sr.size()-1].first;
}

double LookupTable::maximumSth()
{
  return Sth[Sth.size()-1].first;
}

double LookupTable::minimumSr()
{
  return Sr[0].first;
}

double LookupTable::minimumSth()
{
  return Sth[0].first;
}

double LookupTable::maximumLz()
{
  return Lz[Lz.size()-1].first;
}

double LookupTable::minimumLz()
{
  return Lz[0].first;
}

int LookupTable::GrowBlackHole( const double Sr_1, 
                                const double Sth_1, 
                                const double Lz_1, 
                                const double M_1, 
                                const double a_1, 
                                const double M_2, 
                                const double a_2,
                                double consts[4]
                               )
{
  assert(a_2 == this->a);
  assert(M_2 == this->M);
  // First convert Sr_2 and Sth_2 into the same units as our lookup table, namely into units of M_2 = 1; currently M_2 = 500 (for example)
  // The dimensions of Sr and Stheta = [M], since [\sqrt \Theta] = \sqrt{ [M]^2 } and [\int dr] \sqrt{[R]} / [\Delta] = [L] [M]^2/[M]^2 = [M]
  cout << "Min, max, Sr: " << minimumSr() << " " << maximumSr() << " " << Sr_1 << endl;
  cout << "Min, max, Sth: " << minimumSth() << " " << maximumSth() << " " << Sth_1 << endl;
  cout << "Min, max, Lz:  " << minimumLz() << " " << maximumLz() << " " << Lz_1 << endl;

  unordered_map<int,int> keys;
  findKey(Sr_1,  Sr,  keys);
  findKey(Sth_1, Sth, keys);
  int key = findKey(Lz_1,  Lz,  keys);
  if( key == -1 ) return UNSUCCESSFUL; // Did not find any solution
  // Find the corresponding constants of motions
  cout << "Found key with " << key << " and " << keys[key] << endl;
  double Eps_2 = EpsVector[key], Lz_2 =LzVector[key], Q_2 = QVector[key];
  double rmin, rmax;
  assert( abs(helper::I_r(Eps_2,Lz_2,Q_2, a_2,M_2,rmin,rmax) - Sr_1) / abs(Sr_1) < 0.03 );
  double K_2 = Q_2 + pow(Eps_2*a_2 - Lz_2,2);
  // Make sure the limits are OK_2
  double rMin, rMax, thMin, thMax;
  assert( helper::kerr_r_limits (Eps_2,Lz_2,Q_2, a_2, M_2, rMin, rMax));
  assert( helper::kerr_th_limits(Eps_2,Lz_2,Q_2, a_2, M_2, thMin, thMax));
  consts[0] = Eps_2;
  consts[1] = Lz_2;
  consts[2] = Q_2;
  consts[3] = K_2;
  return SUCCESSFUL;
}



bool read_lookup( LookupTable & lookup,
                  const char * filename
                 ) {
  size_t success;
  vector<double> EpsVector;
  vector<double> LzVector;
  vector<double> QVector;
  vector<pair<double,int>> Sr ;
  vector<pair<double,int>> Sth;
  vector<pair<double,int>> Lz;
  int nParticles;
  double rMinLimit, a, M;
  FILE *  f = fopen(filename, "rb");
  // Read number of particles
  success = fread( (char*)&nParticles, sizeof(nParticles), 1, f);
  success = fread( (char*)&rMinLimit,  sizeof(rMinLimit), 1, f);
  success = fread( (char*)&a,  sizeof(a), 1, f);
  success = fread( (char*)&M,  sizeof(M), 1, f);
  assert( nParticles > 100 ); // No way we'll be ok with just 100 entries :)
  assert( rMinLimit > 0 );
  assert( rMinLimit < 1000*M );
  assert( a  <= M && a >= 0 );
  // Resize everything
  EpsVector.resize(nParticles);
  LzVector.resize(nParticles);
  QVector.resize(nParticles);
  Sr.resize(nParticles);
  Sth.resize(nParticles);
  Lz.resize(nParticles);
  assert( EpsVector.size() == nParticles );
  //size_t fread(void *ptr, size_t size_of_elements, size_t number_of_elements, FILE *a_file);
  success = fread( (char*)&EpsVector[0], sizeof(double), EpsVector.size(), f );
  success = fread( (char*)&LzVector[0],  sizeof(double), LzVector.size(),  f );
  success = fread( (char*)&QVector[0],   sizeof(double), QVector.size(),   f );
  // Read Sr, Sth
  {
    vector<double> Sr_1 (Sr.size() );
    vector<int>    Sr_2 (Sr.size() );
    vector<double> Sth_1(Sth.size());
    vector<int>    Sth_2(Sth.size());
    vector<double> Lz_1(Lz.size());
    vector<int>    Lz_2(Lz.size());
    success = fread( (char*)&Sr_1[0],  sizeof(double), Sr_1.size(),  f );
    success = fread( (char*)&Sr_2[0],  sizeof(int),    Sr_2.size(),  f );
    success = fread( (char*)&Sth_1[0], sizeof(double), Sth_1.size(), f );
    success = fread( (char*)&Sth_2[0], sizeof(int),    Sth_2.size(), f );
    success = fread( (char*)&Lz_1[0], sizeof(double), Lz_1.size(), f );
    success = fread( (char*)&Lz_2[0], sizeof(int),    Lz_2.size(), f );
    for( int i =0; i < Sr.size(); ++i ) {
      Sr[i].first  = Sr_1[i] ;
      Sr[i].second = Sr_2[i] ;
      Sth[i].first = Sth_1[i];
      Sth[i].second= Sth_2[i];
      Lz[i].first = Lz_1[i];
      Lz[i].second= Lz_2[i];
    }
  }
  fclose(f);
  // Save into the lookup table
  lookup.EpsVector = EpsVector;
  lookup.LzVector  = LzVector;
  lookup.QVector   = QVector;
  lookup.Sr        = Sr;
  lookup.Sth       = Sth;
  lookup.Lz        = Lz;
  lookup.nParticles= nParticles;
  lookup.rMinLimit = rMinLimit;
  lookup.a         = a;
  lookup.M         = M;
  assert( lookup.EpsVector.size() == EpsVector.size() );
  assert( lookup.Lz.size()        == Lz.size() );
  return true;
}
