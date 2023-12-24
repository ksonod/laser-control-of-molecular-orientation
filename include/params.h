#ifndef params_h
#define params_h

// number of time series
#define SERIES 60000
// maximum number of data
#define NUM 90
// data aquisition -> datum per STEP*dt sec
#define STEP 100

/**molecular properties**/
const double apara = 54.06 * (1.64878*pow(10.0, -41.0)); // α|| in SI unit
const double aperp = 26.09 * (1.64878*pow(10.0, -41.0));// α⊥ in SI unit
const double da = apara - aperp; // Δα
const double bpara = -46.3 * (3.20636 * pow(10.0, -53.0)); //β||
const double bperp = -60.4 * (3.20636 * pow(10.0, -53.0)); //β⊥
const double B = 6081.492475*(3.335641*pow(10.0, -5.0)); // rotational constant in cm^-1
//const double D = 1.301777*pow(10.0,-3.0) * ( 3.335641*pow(10.0,-5.0) ) ; // centrifugal constant in cm^-1
const double D = 0.0; // centrifugal constant in cm^-1
const double Trot = 1.0 / 200.0 / B / vc; // rotational period in sec unit

/*
|    delaym10      |     delay01     |
pulsem1            pulse0            pulse1 (main pulse)
omega              omoega          omega+2omega
*/

const double intensity0 = 20.0*pow(10.0, 12.0); // second pulse (sub pulse) intensity in W/cm^2 unit
const double intensity1 = 30.0*pow(10.0, 12.0); // third pulse (main pulse) intensity in W/cm^2 unit
const double FWHM = 70.0*pow(10.0, -15.0); // FWHM in second unit
const double phase = 0.0;  // relative phase of w and 2w     phi = phase * PI
const double phi = phase * PI; // relative phase
const double rat = 0.5; // intensity ratio I2w/Iw
const double n01 = 0.242; // delay
const double delay01 = n01 * Trot; // delay
int Jmax = 75; //maximum value of J  < NUM


void check_params()
{
    if (Jmax < 3 || Jmax >= NUM)
    {
        throw std::invalid_argument("Jmax should be changed.");
    }
}

#endif /* params_h */
