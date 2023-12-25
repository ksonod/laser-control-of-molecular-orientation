#ifndef runge_kutta_h
#define runge_kutta_h

// Schrodinger equation used for Runge-Kutta calculation.
std::complex<double> calc_schrodinger_equation(double t, int J, int M, std::complex<double> cj3, std::complex<double> cj2, std::complex<double> cj1, std::complex<double> c0, std::complex<double> cJ1, std::complex<double> cJ2, std::complex<double> cJ3){
    double wj3 = (Erot(double(J) - 3.0) - Erot(double(J))) / hbar;
    double wj2 = (Erot(double(J) - 2.0) - Erot(double(J))) / hbar;
    double wj1 = (Erot(double(J) - 1.0) - Erot(double(J))) / hbar;
    double wJ1 = (Erot(double(J) + 1.0) - Erot(double(J))) / hbar;
    double wJ2 = (Erot(double(J) + 2.0) - Erot(double(J))) / hbar;
    double wJ3 = (Erot(double(J) + 3.0) - Erot(double(J))) / hbar;

    double term0 = 0.25 * (E1w(t) * E1w(t) + E2w(t) * E2w(t)) * (aperp + da * d2J0(double(J), double(M)));
    double term1 = (1.0 / 8.0) * (bpara - 3.0*bperp) * E1w(t) * E1w(t) * E2w(t) * d3J3(double(J), double(M)) * cos(phi);
    double term2 = (1.0 / 8.0) * (bpara - 3.0*bperp) * E1w(t) * E1w(t) * E2w(t) * d3j3(double(J), double(M));
    double term3 = 0.25 * da * (E1w(t) * E1w(t) + E2w(t) * E2w(t)) * d2J2(double(J), double(M));
    double term4 = 0.25 * da * (E1w(t) * E1w(t) + E2w(t) * E2w(t)) * d2j2(double(J), double(M));
    double term5 = (3.0 / 8.0 * bperp * E1w(t) * E1w(t) * E2w(t) * d1J1(double(J), double(M)) + (1.0 / 8.0)*(bpara - 3.0*bperp)* E1w(t) * E1w(t) * E2w(t) * d3J1(double(J), double(M))) * cos(phi);
    double term6 = (3.0 / 8.0 * bperp * E1w(t) * E1w(t) * E2w(t) * d1j1(double(J), double(M)) + (1.0 / 8.0)*(bpara - 3.0*bperp)* E1w(t) * E1w(t) * E2w(t) * d3j1(double(J), double(M))) * cos(phi);

    std::complex<double> ret = (I / hbar) * (term0 * c0 + term1 * cJ3 * exp(-I * wJ3*t) + term2 * cj3 * exp(-I * wj3*t) + term3 * cJ2 * exp(-I * wJ2*t) + term4 * cj2 * exp(-I * wj2*t) + term5 * cJ1 * exp(-I * wJ1*t) + term6 * cj1 * exp(-I * wj1*t));
    return  ret;
}


/* This coef function checks if j and m values do no violate rules. Concretely, a physics rule (i.e., J should be larger than |M|) and array restriction (i.e., array index cannot be negative or exceed maximum number of elements) are checked. If everything is fine, an array element is returned.
 */
std::complex<double> coef(std::complex<double> (&c)[num_rot_levels], int j, int m){
    if (j < abs(m)){
        return 0.0;
    }
    else if (j < 0){
        return 0.0;
    }
    else if (num_rot_levels < j){
        return 0.0;
    }
    else{
        return c[j];
    }
}


#endif /* runge_kutta_h */
