/*
LASER-INDUCED MOLECULAR ORIENTATION SIMULATOR
This is a simulator to calculate molecular orientation dynamics induced by two-color intense femtosecond laser pulses.
*/

#include <iostream>
#include <complex>
#include "physics_constant.h"
#include "params.h"
#include "matrix_element.h"
#include "laser_pulses.h"
#include "rotational_energy.h"
#include "runge_kutta.h"


int main()
{
    /**parameters and so on**/
    int Jcalc;// Jcalc is the maximum value of J in calculation
    int N; // number of step
    int step = 1; // step in time evolution
    double cos[SERIES]; // time series of <costheta>
    double cos2[SERIES]; // time series of <cos2theta>
    double NORM = 0.0, Imcos2 = 0.0, Imcos = 0.0; // norm and imaginary part of <cos^2theta> and <costheta> . check these values for the test of the validity of calculation
    double cosmax = 0.0; // maximum value of <costheta>
    double cosmin = 0.0; // minimum value of <costheta>
    double cos2max = 0.0; // maximum value of <cos^2theta>
    double cos20 = 0.0; // <cos^2theta> at t = 0
    clock_t start, end; // calculation time
    double calctime; // calculation time
//    std::complex<double> k[4][NUM], k0[NUM]; //c[] is coefficient of wavefnction.If you change Jmax,you should change size of c[](Jmax/2+1).   k[RK][level]
    std::complex<double> c[NUM], k[4][NUM]; //c[] is coefficient of wavefnction.If you change Jmax,you should change size of c[](Jmax/2+1).
    double cfin[NUM];
    std::complex<double> cj3 = 0.0, cj2 = 0.0, cj1 = 0.0, cJ0 = 0.0, cJ1 = 0.0, cJ2 = 0.0, cJ3 = 0.0;
    std::complex<double> kj3 = 0.0, kj2 = 0.0, kj1 = 0.0, kJ0 = 0.0, kJ1 = 0.0, kJ2 = 0.0, kJ3 = 0.0;

    std::complex<double> align = 0.0, asum = 0.0; // align: alignment parameter
    std::complex<double> orient = 0.0, osum = 0.0; // orient: orientaion parameter

    
    double T = 0.8; // temperature in K unit

    double t;  // time
    double tmin = -2000 * pow(10.0, -15.0) - delay01;
    double tmax = 200.0*pow(10.0, -12.0);
    double dt = 1.0*pow(10.0, -15.0);  // step
    double dtref = dt;
    double dtlong = 40.0*pow(10.0, -15.0);

    double Ethr = 600000.0;
    /*************************************************/
    time_t timer;
    time(&timer); // get the date

    check_params();

    printf("*******************************************************************\n");
    printf(" %s", ctime(&timer));
    printf("*******************************************************************\n");
    printf(" INTENSITY (1st): %f TW \n", intensity0*pow(10.0, -12.0));
    printf(" INTENSITY (2nd): omega1 = %f TW , omega2 = %f TW\n", intensity1*pow(10.0, -12.0), rat*intensity1*pow(10.0, -12.0));
    printf(" PULSE DURATION : FWHM = %f fs\n", FWHM*pow(10.0, 15.0));
    printf(" RELATIVE PHASE : phi = %f*pi\n", phase);

//    FILE *fp1, *fp2, *fp3, *fp4;

    FILE* fp1 = std::fopen("cos2.txt", "w");  // make a text file to write the final results
    FILE* fp2 = std::fopen("cos.txt", "w");
    FILE* fp3 = std::fopen("output.txt", "w");
    FILE* fp4 = std::fopen("pulseenvelope.txt", "w");
    FILE* finpop = std::fopen("finpop.txt", "w");

    // output file
    fprintf(fp3, "*******************************************************************\n");
    fprintf(fp3, " CO_ORIENTATION_VER4-1_131001\n by K.Sonoda\n 01 Oct. 2013\n");
    fprintf(fp3, "*******************************************************************\n");
    fprintf(fp3, "-- DATE --\n");
    fprintf(fp3, " %s", ctime(&timer));
    fprintf(fp3, "*******************************************************************\n");
    fprintf(fp3, "-- VALIDITY --\n");


    /**initial condition**/
    for (N = 0; N < SERIES; N++) // time series of the alignment and orientation parameter
    {
        cos[N] = 0.0;  // set the value of <costheta> zero
        cos2[N] = 0.0;  // set the value of <cos^2theta> zero
    }

    double P; // statistic weight

    Jcalc = calculate_initial_population(T);
    
     /**initialization of final population**/
    for (int J = 0; J<NUM; J++)
    {
        cfin[J] = 0.0;
    }
    // end of the initialization of final population

    printf(" TEMPERATURE    : %f K\n", T);
    printf(" Jcalc          : %d\n", Jcalc);
    printf("*******************************************************************\n");

    /***beginning of the calculation***/
    start = clock(); // calculation time

    for (int Jint = 0; Jint <= Jcalc; Jint++)
    {
        // M calculation
        for (int M = 0; M <= Jint; M++)
        {
            P = calculate_Mlevel_population(T, Jint, M);
            
            N = 0; // set the number of steps 0
            step = 1;

            //initialization of the population
            for (int j = 0; j < NUM; j++)
            {
                //initial population
                if (j == Jint)
                    c[j] = std::complex<double>(1.0, 0.0);
                else
                    c[j] = std::complex<double>(0.0, 0.0);
            }
            //time evolution

            for (t = tmin; t <= tmax; t = t + dt)
            {
                //initialization
                align = 0.0; // align: alignment parameter
                asum = 0.0;
                orient = 0.0; // orient: orientaion parameter
                osum = 0.0;
                NORM = 0.0;
                Imcos2 = 0.0;
                Imcos = 0.0;
                dt = dtref; // small step

                if (Ethr < (E1w(t) + E2w(t))) // Runge-Kutta calculation in the region that the effect of laser pulse is important
                {
                    // calculation of k[0][]
                    for (int j = 0; j <= Jmax; j++)
                    {
//                        cj3 = (0.0, 0.0), cj2 = (0.0, 0.0), cj1 = (0.0, 0.0), cJ0 = (0.0, 0.0), cJ1 = (0.0, 0.0), cJ2 = (0.0, 0.0), cJ3 = (0.0, 0.0);
                        auto [cj3, cj2, cj1, cJ0, cJ1, cJ2, cJ3] = assign_coeff(c, j, M);
                        k[0][j] = dt * RK(t, j, M, cj3, cj2, cj1, cJ0, cJ1, cJ2, cJ3);
                    }// end of the k[0][j] calculation

                     // calculation of k[1][]
                    for (int j = 0; j <= Jmax; j++)
                    {
                        auto[cj3, cj2, cj1, cJ0, cJ1, cJ2, cJ3] = assign_coeff(c, j, M);
                        auto[kj3, kj2, kj1, kJ0, kJ1, kJ2, kJ3] = assign_coeff(k[0], j, M);

                        k[1][j] = dt * RK(t + 0.5 * dt, j, M, cj3 + kj3 * 0.5, cj2 + kj2 * 0.5, cj1 + kj1 * 0.5, cJ0 + kJ0 * 0.5, cJ1 + kJ1 * 0.5, cJ2 + kJ2 * 0.5, cJ3 + kJ3 * 0.5);
                    }// end of the k[1][j] calculation

                     // calculation of k[2][]
                    for (int j = 0; j <= Jmax; j++)
                    {
                        auto[cj3, cj2, cj1, cJ0, cJ1, cJ2, cJ3] = assign_coeff(c, j, M);
                        auto[kj3, kj2, kj1, kJ0, kJ1, kJ2, kJ3] = assign_coeff(k[1], j, M);

                        k[2][j] = dt * RK(t + 0.5 * dt, j, M, cj3 + kj3 * 0.5, cj2 + kj2 * 0.5, cj1 + kj1 * 0.5, cJ0 + kJ0 * 0.5, cJ1 + kJ1 * 0.5, cJ2 + kJ2 * 0.5, cJ3 + kJ3 * 0.5);
                    }// end of the k[2][j] calculation

                     // calculation of k[3][]
                    for (int j = 0; j <= Jmax; j++)
                    {
                        auto[cj3, cj2, cj1, cJ0, cJ1, cJ2, cJ3] = assign_coeff(c, j, M);
                        auto[kj3, kj2, kj1, kJ0, kJ1, kJ2, kJ3] = assign_coeff(k[2], j, M);
                        k[3][j] = dt * RK(t + dt, j, M, cj3 + kj3, cj2 + kj2, cj1 + kj1, cJ0 + kJ0, cJ1 + kJ1, cJ2 + kJ2, cJ3 + kJ3);
                    }// end of the k[3][j] calculation

                    for (int j = 0; j <= Jmax; j++)
                    {
                        c[j] = c[j] + (k[0][j] + 2.0*k[1][j] + 2.0*k[2][j] + k[3][j]) / 6.0;
                    }
                } // end of the RK calculation

                else // long duration step without RK calculation in the region that the laser pulse does not exist
                    dt = dtlong;

                // <cos^2theta> calculation
                for (int j = 0; j <= Jmax; j++)
                {
                    if (j < abs(M)) // impossible
                        asum = 0.0;

                    else if (j == Jmax || j == Jmax - 1)
                    {
                        if (j == abs(M) || j == abs(M) + 1)
                            asum = conj(c[j])*c[j] * d2J0(double(j), double(M));
                        else
                            asum = conj(c[j])*c[j] * d2J0(double(j), double(M)) + conj(c[j])*c[j - 2] * d2j2(double(j), double(M))*exp(-I * (Erot(double(j) - 2.0) - Erot(double(j))) * (t + dt) / hbar);
                        align = asum + align;
                    }

                    else if (j == 0 || j == 1)
                    {
                        asum = conj(c[j])*c[j] * d2J0(double(j), double(M)) + conj(c[j])*c[j + 2] * d2J2(double(j), double(M))*exp(-I * (Erot(double(j) + 2.0) - Erot(double(j))) * (t + dt) / hbar);
                        align = asum + align;
                    }

                    else
                    {
                        if (j - 2 < abs(M))
                            asum = conj(c[j])*c[j] * d2J0(double(j), double(M)) + conj(c[j])*c[j + 2] * d2J2(double(j), double(M))*exp(-I * (Erot(double(j) + 2.0) - Erot(double(j))) * (t + dt) / hbar);
                        else
                            asum = conj(c[j])*c[j] * d2J0(double(j), double(M)) + conj(c[j])*c[j - 2] * d2j2(double(j), double(M))*exp(-I * (Erot(double(j) - 2.0) - Erot(double(j))) * (t + dt) / hbar) + conj(c[j])*c[j + 2] * d2J2(double(j), double(M))*exp(-I * (Erot(double(j) + 2.0) - Erot(double(j))) * (t + dt) / hbar);
                        align = asum + align;
                    }
                }// end of the <cos^2theta> calculation

                 //imaginary part of <cos^2theta>
                Imcos2 = imag(align);

                // <costheta> calculation
                for (int j = 0; j <= Jmax; j++)
                {
                    if (j < abs(M))
                        osum = 0.0;

                    else if (j == Jmax) // no interaction with upper state
                    {
                        if (j == abs(M)) // no interaction
                            osum = 0.0;
                        else // interaction with lower state(J=J-1)
                            osum = conj(c[j])*c[j - 1] * d1j1(double(j), double(M))*exp(-I * (Erot(double(j) - 1.0) - Erot(double(j))) * (t + dt) / hbar);
                        orient = osum + orient;
                    }

                    else if (j == 0) // no interaction with lower state
                    {
                        osum = conj(c[j])*c[j + 1] * d1J1(double(j), double(M))*exp(-I * (Erot(double(j) + 1.0) - Erot(double(j))) * (t + dt) / hbar);
                        orient = osum + orient;
                    }

                    else
                    {
                        if (j == abs(M)) // no interaction with lower state
                            osum = conj(c[j])*c[j + 1] * d1J1(double(j), double(M))*exp(-I * (Erot(double(j) + 1.0) - Erot(double(j))) * (t + dt) / hbar);
                        else
                            osum = conj(c[j])*c[j - 1] * d1j1(double(j), double(M))*exp(-I * (Erot(double(j) - 1.0) - Erot(double(j))) * (t + dt) / hbar) + conj(c[j])*c[j + 1] * d1J1(double(j), double(M))*exp(-I * (Erot(double(j) + 1.0) - Erot(double(j))) * (t + dt) / hbar);
                        orient = osum + orient;
                    }
                }// end of the <costheta> calculation

                 //imaginary part of <costheta>
                Imcos = imag(orient);

                //norm calculation
                for (int j = 0; j <= Jmax; j++)
                {
                    NORM = norm(c[j]) + NORM;
                }// end of the norm calculation

                if (T == 0.0)
                {
                    printf("%f\t%f\t%f\t%f\t%f\n", (t + dt)*pow(10.0, 15.0), real(orient), abs(Imcos), abs(Imcos2), NORM);
                }

                // get alignment and orientation parameters
                if (Ethr < (E1w(t) + E2w(t))) // at small step RK calculation
                {
                    fprintf(fp4, "%f\t%f\n", (t + dt) / Trot, 0.5* (E1w(t + dt) + E2w(t + dt))*(E1w(t + dt) + E2w(t + dt))*epsilon*vc * pow(10.0, -12.0)* pow(10.0, -4.0));

                    if (step%STEP == 0 && STEP <= step) // data aquisition per #STEP steps
                    {
                        cos2[N] = P * real(align) + cos2[N];  // time series of <cos^2theta>
                        cos[N] = P * real(orient) + cos[N];  // time series of <costheta>

                        // at the end of the calculation...
                        if (Jint == Jcalc && M == Jint)
                        {
                            // write the results in txt file
                            //                            fprintf(fp1,"%f\t%f\n",(t+dt)*pow(10.0,15.0),cos2[N]) ;
                            //                            fprintf(fp2,"%f\t%f\n",(t+dt)*pow(10.0,15.0),cos[N]) ;
                            fprintf(fp1, "%f\t%f\n", (t + dt) / Trot, cos2[N]);
                            fprintf(fp2, "%f\t%f\n", (t + dt) / Trot, cos[N]);
                            //                            fprintf(fp4, "%f\t%f\n", (t + dt) / Trot, E1w(t+dt)+ E2w(t + dt));

                            // obtain the maximum value of <cos^2theta>
                            if (cos2max <= cos2[N])
                                cos2max = cos2[N];

                            // obtain the value of <cos^2theta> around t = 0
                            if ((-STEP * dtref / 2.0 < t) && (t < STEP*dtref / 2.0))
                                cos20 = cos2[N];

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

                    else;

                    step = step + 1;

                }

                else // at long duration step calculation
                {
                    cos2[N] = P * real(align) + cos2[N];  // time series of <cos^2theta>
                    cos[N] = P * real(orient) + cos[N];  // time series of <costheta>

                    fprintf(fp4, "%f\t%f\n", (t + dt) / Trot, 0.5* (E1w(t + dt) + E2w(t + dt)) * (E1w(t + dt) + E2w(t + dt)) * epsilon * vc * pow(10.0, -12.0) * pow(10.0, -4.0));

                    // at the end of the calculation...
                    if (Jint == Jcalc && M == Jint)
                    {
                        // write the results in txt file
                        //                        fprintf(fp1,"%f\t%f\n",(t+dt)*pow(10.0,15.0),cos2[N]) ;
                        //                        fprintf(fp2,"%f\t%f\n",(t+dt)*pow(10.0,15.0),cos[N]) ;
                        fprintf(fp1, "%f\t%f\n", (t + dt) / Trot, cos2[N]);
                        fprintf(fp2, "%f\t%f\n", (t + dt) / Trot, cos[N]);
                        //                        fprintf(fp4, "%f\t%f\n", (t + dt) / Trot, E1w(t + dt) + E2w(t + dt));

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

                // avoid fatal error
                if (SERIES <= N) // if N exceed the number of SERIES, kill the program.
                {
                    printf("\n\n - ERROR! - \n\n");
                    return 1;
                }

            }// end of the time evolution

             // final population calculation
            for (int j = 0; j <= Jmax; j++)
            {
                cfin[j] = P * norm(c[j]) + cfin[j];
                if (Jint == Jcalc && M == Jint) // at the end of the calculation, write the results in txt file
                    fprintf(finpop, "%d\t%f\n", j, cfin[j]);
            }// end of the final population calculation

            printf("J%dM%d\tN=%f Im(cos)=%f Im(cos2)=%f\n", Jint, M, NORM, abs(Imcos), abs(Imcos2));
            fprintf(fp3, " J%dM%d\tN=%f Im(cos)=%f Im(cos2)=%f\n", Jint, M, NORM, abs(Imcos), abs(Imcos2));

        }// end of the calculation of all M for specific J

    } // end of the all Jint calculation

    end = clock(); // calculation time
    calctime = (double(end - start)) / CLOCKS_PER_SEC; // calculation time (seconds)

    // output the maximum value of <costheta>
    printf("-----------------------\n");
    printf("calculation time = %f sec\n", calctime);
    printf("<cos>max  = %f\n", cosmax);
    printf("<cos>min  = %f\n", cosmin);
    printf("<cos2>max = %f\n", cos2max);
    printf("<cos2>t=0 = %f\n", cos20);


    // write the conditions and results in txt file
    fprintf(fp3, "*******************************************************************\n");
    fprintf(fp3, "-- LASER PULSE --\n");
    fprintf(fp3, " -first pulse\n");
    fprintf(fp3, " INTENSITY      : %f TW\n", intensity0*pow(10.0, -12.0));
    fprintf(fp3, " PULSE DURATION : FWHM = %f fs\n", FWHM*pow(10.0, 15.0));
    //    fprintf(fp3," DELAY          : %f ps (%f*Trot)\n", delay*pow(10.0,12.0) , n );
    fprintf(fp3, " -second pulse\n");
    fprintf(fp3, " INTENSITY      : omega1 = %f TW , omega2 = %f TW\n", intensity1*pow(10.0, -12.0), rat*intensity1*pow(10.0, -12.0));
    fprintf(fp3, " PULSE DURATION : FWHM = %f fs\n", FWHM*pow(10.0, 15.0));
    fprintf(fp3, " RELATIVE PHASE : phi = %f*pi\n", phase);
    fprintf(fp3, "*******************************************************************\n");
    fprintf(fp3, "-- TEMPERATURE --\n");
    fprintf(fp3, " TEMPERATURE    : %f K\n", T);
    fprintf(fp3, "*******************************************************************\n");
    fprintf(fp3, "-- ROTATIONAL LEVEL --\n");
    fprintf(fp3, " Jcalc          : %d\n", Jcalc);
    fprintf(fp3, " Jmax           : %d\n", Jmax);
    fprintf(fp3, "*******************************************************************\n");
    fprintf(fp3, "-- TIME --\n");
    fprintf(fp3, " tmin           : %f fs\n", tmin*pow(10.0, 15.0));
    fprintf(fp3, " tmax           : %f fs\n", tmax*pow(10.0, 15.0));
    fprintf(fp3, " dt             : %f fs\n", dtref*pow(10.0, 15.0));
    fprintf(fp3, " dtlong         : %f fs\n", dtlong*pow(10.0, 15.0));
    fprintf(fp3, "*******************************************************************\n");
    fprintf(fp3, "-- RESULTS --\n");
    fprintf(fp3, " <cos>max       : %f\n", cosmax);
    fprintf(fp3, " <cos>min       : %f\n", cosmin);
    fprintf(fp3, " <cos^2>max     : %f\n", cos2max);
    fprintf(fp3, " <cos^2>t=0     : %f\n", cos20);
    fprintf(fp3, "*******************************************************************\n");
    fprintf(fp3, "-- CALCULATION TIME --\n");
    fprintf(fp3, " CALC. TIME     : %f sec\n", calctime);
    fprintf(fp3, "*******************************************************************\n");
    
    std::fclose(fp1);
    std::fclose(fp2);
    std::fclose(fp3);
    std::fclose(fp4);
    std::fclose(finpop);

    return 0;

}
