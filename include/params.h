#ifndef params_h
#define params_h

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
double T = 0.8; // temperature in K unit


const double intensity0 = 20.0*pow(10.0, 12.0); // second pulse (sub pulse) intensity in W/cm^2 unit
const double intensity1 = 30.0*pow(10.0, 12.0); // third pulse (main pulse) intensity in W/cm^2 unit
const double FWHM = 70.0*pow(10.0, -15.0); // FWHM in second unit
const double phase = 0.0;  // relative phase of w and 2w     phi = phase * PI
const double phi = phase * PI; // relative phase
const double rat = 0.5; // intensity ratio I2w/Iw
const double n_delay = 0.242; // delay
const double t_delay = n_delay * Trot; // delay
const int Jmax = 75; //maximum value of J  < num_rot_levels
const int num_rot_levels = 90;

double tmin = -2000 * pow(10.0, -15.0) - t_delay;
double tmax = 200.0*pow(10.0, -12.0);
double dt_small = 1.0*pow(10.0, -15.0);  // step
double dt_large = 40.0*pow(10.0, -15.0);
double Ethr = 600000.0;

const int num_time_series_data = 60000;
const int data_sampling_step = 50; // Data is sampled every data_sampling_step*dt seconds.


void check_params()
{
    if (Jmax < 3 || Jmax >= num_rot_levels)
    {
        throw std::invalid_argument("Jmax should be changed.");
    }
}

#endif /* params_h */
