#ifndef runge_kutta_h
#define runge_kutta_h


std::complex<double> RK(double t, int J, int M, std::complex<double> Cj3, std::complex<double> Cj2, std::complex<double> Cj1, std::complex<double> C0, std::complex<double> CJ1, std::complex<double> CJ2, std::complex<double> CJ3)
{
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

    std::complex<double> ret = (I / hbar) * (term0 * C0 + term1 * CJ3 * exp(-I * wJ3*t) + term2 * Cj3 * exp(-I * wj3*t) + term3 * CJ2 * exp(-I * wJ2*t) + term4 * Cj2 * exp(-I * wj2*t) + term5 * CJ1 * exp(-I * wJ1*t) + term6 * Cj1 * exp(-I * wj1*t));
    return  ret;
}


/*
 This funcition considers allowed interactions between different J states.
 */
std::tuple<std::complex<double>, std::complex<double>, std::complex<double>, std::complex<double>, std::complex<double>, std::complex<double>, std::complex<double>> assign_coeff(std::complex<double> (&c)[NUM], int j, int M)
{
    std::complex<double> cj3 = (0.0, 0.0), cj2 = (0.0, 0.0), cj1 = (0.0, 0.0), cJ0 = (0.0, 0.0), cJ1 = (0.0, 0.0), cJ2 = (0.0, 0.0), cJ3 = (0.0, 0.0);

    if (j >= abs(M))
    {
        cJ0 = c[j];

        if (2 < j)
        {
            if (j == (Jmax - 1))
            {
                cJ1 = c[j+1];
            }
            else if (j == (Jmax - 2))
            {
                cJ1 = c[j+1];
                cJ2 = c[j+2];
            }
            else;
            {
                cJ1 = c[j+1];
                cJ2 = c[j+2];
                cJ3 = c[j+3];
            }


            if (j == (abs(M) + 1))
            {
                cj1 = c[j-1];
            }
            else if (j == (abs(M) + 2))
            {
                cj1 = c[j-1];
                cj2 = c[j-2];
            }
            else;
            {
                cj1 = c[j-1];
                cj2 = c[j-2];
                cj3 = c[j-3];
            }
        }
        else  // j = 0, 1, 2
        {
            if (j == 2)
                if (abs(M) == 2)
                {
                    cJ1 = c[j+1];
                    cJ2 = c[j+2];
                    cJ3 = c[j+3];
                }
                else if (abs(M) == 1)
                {
                    cj1 = c[j-1];
                    cJ1 = c[j+1];
                    cJ2 = c[j+2];
                    cJ3 = c[j+3];
                }
                else if (M == 0)
                {
                    cj1 = c[j-1];
                    cj2 = c[j-2];
                    cJ1 = c[j+1];
                    cJ2 = c[j+2];
                    cJ3 = c[j+3];
                }
                else;
            else if (j == 1)
                if (abs(M) == 1)
                {
                    cJ1 = c[j+1];
                    cJ2 = c[j+2];
                    cJ3 = c[j+3];
                }
                else if (M == 0)
                {
                    cj1 = c[j-1];
                    cJ1 = c[j+1];
                    cJ2 = c[j+2];
                    cJ3 = c[j+3];
                }
                else;
            else  // j == 0
            {
                cJ1 = c[j+1];
                cJ2 = c[j+2];
                cJ3 = c[j+3];
            }
        }
    }

    return {cj3, cj2, cj1, cJ0, cJ1, cJ2, cJ3};
}

#endif /* runge_kutta_h */
