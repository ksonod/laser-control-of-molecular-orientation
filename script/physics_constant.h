#include <complex>
#include <math.h>
using namespace std;

#ifndef physics_constant_h
#define physics_constant_h

const complex<double>I = complex<double>(0.0, 1.0);  // imaginary unit
const double PI = acos(-1.0);  // pi
const double kB = 1.3806504*pow(10.0, -23.0); // Boltzmann constant
const double h = 6.62606896*pow(10.0, -34.0); // Planck constant in J•sec
const double hbar = h / 2.0 / PI; // Dirac constant difined as hbar=h/pi
const double vc = 299792458.0; // velocity of light in m/sec
const double epsilon = 8.854187817*pow(10.0, -12.0); // permittivity of vacuum



#endif /* physics_constant_h */
