#ifndef physics_constant_h
#define physics_constant_h

const std::complex<double>I = std::complex<double>(0.0, 1.0);  // imaginary unit
const double pi = acos(-1.0);  // pi
const double k_boltzmann = 1.3806504 * pow(10.0, -23.0); // Boltzmann constant
const double h_planck = 6.62606896 * pow(10.0, -34.0); // Planck constant in J•sec
const double hbar = h_planck / 2.0 / pi; // Dirac constant difined as hbar=h/pi
const double speed_of_light = 299792458.0; // velocity of light in m/sec
const double epsilon = 8.854187817 * pow(10.0, -12.0); // permittivity of vacuum

#endif /* physics_constant_h */
