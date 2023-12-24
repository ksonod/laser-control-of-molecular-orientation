#ifndef laser_pulses_h
#define laser_pulses_h


double pulse_envelope(double t, double intensity, double fwhm)
{
    double envelope = exp(-2.0 * log(2.0) * t * t / fwhm / fwhm);
    double amplitude = sqrt(2.0 * 10000.0 * intensity / epsilon / vc);
    return envelope * amplitude;
}


double E1w(double t) //laser pulse
{
                                // pulse0 -> sub pulse (second pulse)
    double pulse0 = pulse_envelope(t + t_delay, intensity0, FWHM);
//    double envelope0 = exp(-2.0*log(2.0)*(t + delay01)*(t + delay01) / FWHM / FWHM);
//    double E0 = sqrt(2.0 * 10000.0 * intensity0 / epsilon / vc); // amplitude of 800 nm fundumental
//
//    pulse0 = E0 * envelope0; //  t_pulse0 < t_pulse1

                             // pulse1 -> main pulse (second pulse) overlaped with the second harmonic.
    double pulse1 = pulse_envelope(t, intensity1, FWHM);

    return pulse0 + pulse1;
}


double E2w(double t) //laser pulse  (2omega)
{
    return pulse_envelope(t, intensity1*rat, FWHM);

//    double INTENSITYw2 = rat * intensity1; // intensity of SHG
//    double E2 = sqrt(2.0 * 10000.0 * INTENSITYw2 / epsilon / vc); // amplitude of SHG

}

#endif /* laser_pulses_h */
