#ifndef params_h
#define params_h

/* Molecular properties */
const double apara = 54.06 * (1.64878 * pow(10.0, -41.0));  // Parallel component of polarizability in SI unit
const double aperp = 26.09 * (1.64878 * pow(10.0, -41.0));  // Perpendicular component of polarizability in SI unit
const double da = apara - aperp;  // Difference of parallel and perpendicular component of polarizability
const double bpara = -46.3 * (3.20636 * pow(10.0, -53.0));  // Parallel component of hyperpolarizability in SI unit
const double bperp = -60.4 * (3.20636 * pow(10.0, -53.0));  // Perpendicular component of hyperpolarizability in SI unit
const double B_rot = 6081.492475 * (3.335641 * pow(10.0, -5.0));  // Rotational constant in cm^-1
const double D_centr = 0.0;  // Centrifugal distortion constant in cm^-1
const double t_rot_period = 1.0 / 200.0 / B_rot / speed_of_light;  // Rotational period in seconds
double T = 0.8;  // Temperature in K unit

/* Laser parameters */
const double intensity0 = 20.0 * pow(10.0, 12.0); // First laser pulse intensity in W/cm^2 unit
const double intensity1 = 30.0 * pow(10.0, 12.0); // Second laser pulse (fundamental frequency) intensity in W/cm^2 unit
const double pulse_fwhm = 70.0 * pow(10.0, -15.0);  // Laser pulse width (full-width at half maximum) in seconds
const double phase = 0.0;  // Relative phase of w (fundamental) and 2w (second harmonic). ->  phi = phase * PI
const double phi = phase * PI; // Relative phase in radian
const double rat = 0.5;  // Intensity ratio of second harmonic and fundamental light: I2w/Iw
const double n_delay = 0.242; // Temporal delay between first and second pulses in the unit of t_rot_period (rotational period)
const double t_delay = n_delay * t_rot_period;  // Temporal delay between the first and second pulses in seconds.

/* Simulation parameters */
const double tmin = -2000 * pow(10.0, -15.0) - t_delay;  // Minimum time
const double tmax = 200.0 * pow(10.0, -12.0);  // Maximum time
const double dt_small = 1.0 * pow(10.0, -15.0);  // Small time step for the Runge-Kutta method.
const double dt_large = 40.0 * pow(10.0, -15.0);  // Large time step for the region where laser pulses are not strong.
const int Jmax = 75;  // Maximum value of the rotational quantum number, J. Jmax < num_rot_levels
const int num_rot_levels = 90;  // Number of rotational levels, which is needed for creating arrays.
const int num_time_series_data = 60000;  // Number of time-series data points
const int data_sampling_step = 50;  // Data is sampled every data_sampling_step*dt seconds.
const double rot_population_thr = 0.0001;   // threshold of rotationl distribution
const double electric_field_thr = 600000.0;  // Threshold value for electric field strength to judge if the light is strong. -> Decision point for RK calculation.


void check_params(){
    if (Jmax < 3 || num_rot_levels <= Jmax){
        throw std::invalid_argument("Jmax should be changed.");
    }
}

#endif /* params_h */
