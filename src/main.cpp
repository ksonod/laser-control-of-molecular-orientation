/*
LASER-INDUCED MOLECULAR ORIENTATION SIMULATOR
This is a simulator to calculate molecular orientation dynamics induced by two-color intense femtosecond laser pulses.
*/

#include <iostream>
#include <fstream>
#include <complex>
#include "physics_constant.h"
#include "params.h"
#include "matrix_element.h"
#include "laser_pulses.h"
#include "rotational_energy.h"
#include "runge_kutta.h"
#include "expectation_values.h"


int main()
{
    /**parameters and so on**/
    double P; // statistic weight
    int Jcalc;// Jcalc is the maximum value of J in calculation
    int N; // number of step
    int step = 1; // step in time evolution
    double cos[SERIES] = {0}; // time series of <costheta>
    double cos2[SERIES] = {0}; // time series of <cos2theta>
    double NORM = 0.0, Imcos2 = 0.0, Imcos = 0.0; // norm and imaginary part of <cos^2theta> and <costheta> . check these values for the test of the validity of calculation
    double cosmax = 0.0; // maximum value of <costheta>
    double cosmin = 0.0; // minimum value of <costheta>
    double cos2max = 0.0; // maximum value of <cos^2theta>
    double cos20 = 0.0; // <cos^2theta> at t = 0
    clock_t start, end; // calculation time
    double calctime; // calculation time
//    std::complex<double> k[4][NUM], k0[NUM]; //c[] is coefficient of wavefnction.If you change Jmax,you should change size of c[](Jmax/2+1).   k[RK][level]
    std::complex<double> c[NUM], k[4][NUM]; //c[] is coefficient of wavefnction.If you change Jmax,you should change size of c[](Jmax/2+1).
    double cfin[NUM] = {0};
    std::complex<double> cj3, cj2, cj1, cJ0, cJ1, cJ2, cJ3;
    std::complex<double> kj3, kj2, kj1, kJ0, kJ1, kJ2, kJ3;

    std::complex<double> align = 0.0, asum = 0.0; // align: alignment parameter
    std::complex<double> orient = 0.0, osum = 0.0; // orient: orientaion parameter
    
    double T = 0.8; // temperature in K unit

    double t;  // time
    double dt = dt_small;
    /*************************************************/
    time_t timer;
    time(&timer); // get the date

    check_params();
    Jcalc = calculate_initial_population(T);

    std::cout << "*******************************************************************\n";
    std::cout << ctime(&timer);
    std::cout << "*******************************************************************\n";
    std::cout << " INTENSITY (1st): " << intensity0*pow(10.0, -12.0) << " TW \n";
    std::cout << " INTENSITY (2nd): omega1 = " << intensity1*pow(10.0, -12.0) << " TW , omega2 = " << rat*intensity1*pow(10.0, -12.0) << " TW\n";
    std::cout << " PULSE DURATION : FWHM = " << FWHM*pow(10.0, 15.0) << " fs\n";
    std::cout << " RELATIVE PHASE : phi = " << phase << "*pi\n";
    std::cout << " TEMPERATURE    : " << T << " K\n";
    std::cout << " Jcalc          : " << Jcalc << "\n";
    std::cout << "*******************************************************************\n";

    std::ofstream file_cos2("cos2.txt");
    std::ofstream file_cos("cos.txt");
    std::ofstream file_output("output.txt");
    std::ofstream file_finpop("finpop.txt");
    
    
    // output file
    file_output << "*******************************************************************\n";
    file_output << "-- DATE --\n";
    file_output << ctime(&timer) << "\n";
    file_output << "*******************************************************************\n";
    file_output << "-- LOG --\n";

    /*beginning of the calculation*/
    start = clock(); // calculation time

    for (int Jint = 0; Jint <= Jcalc; Jint++){
        // M calculation
        for (int M = 0; M <= Jint; M++){
            P = calculate_Mlevel_population(T, Jint, M);
            
            N = 0; // set the number of steps 0
            step = 1;

            //initialization of the population
            for (int j = 0; j < NUM; j++){
                //initial population
                if (j == Jint)
                    c[j] = std::complex<double>(1.0, 0.0);
                else
                    c[j] = std::complex<double>(0.0, 0.0);
            }
            //time evolution

            for (t = tmin; t <= tmax; t = t + dt){
                //initialization
                align = 0.0; // alignment parameter
                asum = 0.0;
                orient = 0.0; // orientaion parameter
                osum = 0.0;
                NORM = 0.0;
                Imcos2 = 0.0;
                Imcos = 0.0;
                dt = dt_small; // small step

                if (Ethr < (E1w(t) + E2w(t))){  // Runge-Kutta calculation in the region that the effect of laser pulse is important
                    // calculation of k[0][]
                    for (int j = 0; j <= Jmax; j++){
                        auto [cj3, cj2, cj1, cJ0, cJ1, cJ2, cJ3] = assign_coeff(c, j, M);
                        k[0][j] = dt * RK(t, j, M, cj3, cj2, cj1, cJ0, cJ1, cJ2, cJ3);
                    }

                     // calculation of k[1][]
                    for (int j = 0; j <= Jmax; j++){
                        auto[cj3, cj2, cj1, cJ0, cJ1, cJ2, cJ3] = assign_coeff(c, j, M);
                        auto[kj3, kj2, kj1, kJ0, kJ1, kJ2, kJ3] = assign_coeff(k[0], j, M);

                        k[1][j] = dt * RK(t + 0.5 * dt, j, M, cj3 + kj3 * 0.5, cj2 + kj2 * 0.5, cj1 + kj1 * 0.5, cJ0 + kJ0 * 0.5, cJ1 + kJ1 * 0.5, cJ2 + kJ2 * 0.5, cJ3 + kJ3 * 0.5);
                    }

                     // calculation of k[2][]
                    for (int j = 0; j <= Jmax; j++){
                        auto[cj3, cj2, cj1, cJ0, cJ1, cJ2, cJ3] = assign_coeff(c, j, M);
                        auto[kj3, kj2, kj1, kJ0, kJ1, kJ2, kJ3] = assign_coeff(k[1], j, M);

                        k[2][j] = dt * RK(t + 0.5 * dt, j, M, cj3 + kj3 * 0.5, cj2 + kj2 * 0.5, cj1 + kj1 * 0.5, cJ0 + kJ0 * 0.5, cJ1 + kJ1 * 0.5, cJ2 + kJ2 * 0.5, cJ3 + kJ3 * 0.5);
                    }

                     // calculation of k[3][]
                    for (int j = 0; j <= Jmax; j++){
                        auto[cj3, cj2, cj1, cJ0, cJ1, cJ2, cJ3] = assign_coeff(c, j, M);
                        auto[kj3, kj2, kj1, kJ0, kJ1, kJ2, kJ3] = assign_coeff(k[2], j, M);
                        k[3][j] = dt * RK(t + dt, j, M, cj3 + kj3, cj2 + kj2, cj1 + kj1, cJ0 + kJ0, cJ1 + kJ1, cJ2 + kJ2, cJ3 + kJ3);
                    }

                    for (int j = 0; j <= Jmax; j++){
                        c[j] = c[j] + (k[0][j] + 2.0*k[1][j] + 2.0*k[2][j] + k[3][j]) / 6.0;
                        NORM += norm(c[j]);
                    }
                } // end of the RK calculation
                else{ // Large step without RK calculation in the region where the laser pulse is extremely weak or does not exist.
                
                    dt = dt_large;
                
                    // Norm calculation
                    for (int j = 0; j <= Jmax; j++)
                    {
                        NORM += norm(c[j]);
                    }
                }

                align = calculate_cos2_expectation_value(c, M, t, dt);
                orient = calculate_cos_expectation_value(c, M, t, dt);

                if (T == 0.0)
                {
                    std::cout << (t + dt)*pow(10.0, 15.0) << "\t" << real(orient) << "\t" << orient.imag() << "\t" << align.imag() << "\t" << NORM << "\n";
                }

                // get alignment and orientation parameters
                if (Ethr < (E1w(t) + E2w(t))) // at small step RK calculation
                {
                    if (step%STEP == 0 && STEP <= step) // data aquisition per #STEP steps
                    {
                        cos2[N] = P * real(align) + cos2[N];  // time series of <cos^2theta>
                        cos[N] = P * real(orient) + cos[N];  // time series of <costheta>

                        // at the end of the calculation...
                        if (Jint == Jcalc && M == Jint)
                        {
                            // write the results in txt file
                            file_cos2 << (t + dt) / Trot << "\t" << cos2[N] << "\n";
                            file_cos << (t + dt) / Trot << "\t" << cos[N] << "\n";

                            // obtain the maximum value of <cos^2theta>
                            if (cos2max <= cos2[N])
                                cos2max = cos2[N];

                            // obtain the value of <cos^2theta> around t = 0
                            if ((-STEP * dt_small / 2.0 < t) && (t < STEP*dt_small / 2.0))
                                cos20 = cos2[N];

                            // obtain the maximum value of <costheta>
                            if (cosmax <= cos[N])
                                cosmax = cos[N];

                            // obtain the minimum value of <costheta>
                            if (cos[N] <= cosmin)
                                cosmin = cos[N];
                        }

                        else;

                        N += 1; //next step data
                    }

                    else;

                    step += 1;

                }

                else // at long duration step calculation
                {
                    cos2[N] = P * real(align) + cos2[N];  // time series of <cos^2theta>
                    cos[N] = P * real(orient) + cos[N];  // time series of <costheta>

                    // at the end of the calculation...
                    if (Jint == Jcalc && M == Jint)
                    {
                        // write the results in txt file
                        file_cos2 << (t + dt) / Trot << "\t" << cos2[N] << "\n";
                        file_cos << (t + dt) / Trot << "\t" << cos[N] << "\n";

                        // obtain the maximum value of <cos^2theta>
                        if ((0 <= t) && (cos2max <= cos2[N]))
                            cos2max = cos2[N];

                        // obtain the maximum value of <costheta>
                        if (cosmax <= cos[N])
                            cosmax = cos[N];

                        // obtain the minimum value of <costheta>
                        if (cos[N] <= cosmin)
                            cosmin = cos[N];
                    }
                    else;

                    N = N + 1; //next step data
                }
                // end of the getting the alignment and orientation parameters
            }// end of the time evolution

             // final population calculation
            for (int j = 0; j <= Jmax; j++)
            {
                cfin[j] = P * norm(c[j]) + cfin[j];
                if (Jint == Jcalc && M == Jint) // at the end of the calculation, write the results in txt file
                    file_finpop << j << "\t" << cfin[j] << "\n";
            }// end of the final population calculation

            std::cout << "J" << Jint << "M" << M << "\tN=" << NORM << "\tIm(cos)=" << orient.imag() << "\tIm(cos2)=" << align.imag() << "\n";
            file_output << " J" << Jint << "M" << M << "\tN=" << NORM << "\tIm(cos)=" << orient.imag() << "\tIm(cos2)=" << orient.imag() << "\n";

        }// end of the calculation of all M for specific J
        
    } // end of the all Jint calculation

    end = clock(); // calculation time
    calctime = (double(end - start)) / CLOCKS_PER_SEC; // calculation time (seconds)

    // output the maximum value of <costheta>
    std::cout << "-----------------------\n";
    std::cout << "calculation time = " << calctime << " sec.\n";
    std::cout << "<cos>max  = " << cosmax << "\n";
    std::cout << "<cos>min  = " << cosmin << "\n";
    std::cout << "<cos2>max = " << cos2max << "\n";
    std::cout << "<cos2>t=0 = " << cos20 << "\n";

    // write the conditions and results in txt file
    file_output << "*******************************************************************\n";
    file_output << "-- LASER PULSE --\n";
    file_output << " -first pulse\n";
    file_output << " INTENSITY      : " << intensity0*pow(10.0, -12.0) << " TW\n";
    file_output << " PULSE DURATION : FWHM = " << FWHM*pow(10.0, 15.0) << " fs\n";
    file_output << " -second pulse\n";
    file_output << " INTENSITY      : omega1 = " << intensity1*pow(10.0, -12.0) << " TW , omega2 = " << rat*intensity1*pow(10.0, -12.0) << " TW\n";
    file_output << " PULSE DURATION : FWHM = " << FWHM*pow(10.0, 15.0) << " fs\n";
    file_output << " RELATIVE PHASE : phi = " << phase << "*pi\n";
    file_output << " DELAY          : " <<  t_delay*pow(10.0,12.0) << " ps (" << n_delay << "Trot)\n";
    file_output << "*******************************************************************\n";
    file_output << "-- TEMPERATURE --\n";
    file_output << " TEMPERATURE    : " << T << " K\n";
    file_output << "*******************************************************************\n";
    file_output << "-- ROTATIONAL LEVEL --\n";
    file_output << " Jcalc          : " << Jcalc << "\n";
    file_output << " Jmax           : " << Jmax << "\n";
    file_output << "*******************************************************************\n";
    file_output << "-- TIME --\n";
    file_output << " tmin           : " << tmin*pow(10.0, 15.0) << " fs\n";
    file_output << " tmax           : " << tmax*pow(10.0, 15.0) << " fs\n";
    file_output << " dt             : " << dt_small*pow(10.0, 15.0) << " fs\n";
    file_output << " dtlong         : " << dt_large*pow(10.0, 15.0) << " fs\n";
    file_output << "*******************************************************************\n";
    file_output << "-- RESULTS --\n";
    file_output << " <cos>max       : " << cosmax << " \n";
    file_output << " <cos>min       : " << cosmin << " \n";
    file_output <<  " <cos^2>max     : " << cos2max << "\n";
    file_output << "*******************************************************************\n";
    file_output << "-- CALCULATION TIME --\n";
    file_output << " CALC. TIME     : " << calctime << " sec\n";
    file_output << "*******************************************************************\n";
    
    file_cos2.close();
    file_cos.close();
    file_finpop.close();

    return 0;

}
