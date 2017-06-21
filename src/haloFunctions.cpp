#include <cmath>
double dM(int n, double a, double b) {
  return (4.*(pow(a,3 - n) - pow(b,3 - n))*M_PI)/(-3. + n);
}