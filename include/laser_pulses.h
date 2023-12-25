#ifndef laser_pulses_h
#define laser_pulses_h

// Gaussian pulse envelope
double pulse_envelope(double t, double intensity, double fwhm){
    double envelope = exp(-2.0 * log(2.0) * t * t / fwhm / fwhm);
    double amplitude = sqrt(2.0 * 10000.0 * intensity / epsilon / speed_of_light);
    return envelope * amplitude;
}


// Laser pulse with the fundamental frequency
double E1w(double t){
    // First laser pulse
    double pulse0 = pulse_envelope(t + t_delay, intensity0, pulse_fwhm);

    // Second pulse overlaped with the second harmonic
    double pulse1 = pulse_envelope(t, intensity1, pulse_fwhm);

    return pulse0 + pulse1;
}


// Laser pulse with the second harmonic
double E2w(double t){
    return pulse_envelope(t, intensity1 * rat, pulse_fwhm);
}

#endif /* laser_pulses_h */
