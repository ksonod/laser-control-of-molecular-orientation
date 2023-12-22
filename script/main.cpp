/*
Molecular Orientation Simulator
 
This is a simulator to calculate laser-induced molecular orientation dynamics.
*/


#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <complex>
#include <math.h>
#include <time.h>
#include "physics_constant.h"
#include "matrix_element.h"
#include "laser_pulses.h"

using namespace std;


// number of time series
#define SERIES 60000
// maximum number of data
#define NUM 90
// data aquisition -> datum per STEP*dt sec
#define STEP 100

                                     /**functions**/
double d3J3(double J, double M); // <cos^3> delta J = +3
double d3J1(double J, double M); // <cos^3> delta J = +1
double d3j1(double J, double M); // <cos^3> delta J = -1
double d3j3(double J, double M); // <cos^3> delta J = -3
double d2J2(double J, double M); // <cos^2> delta J = +2
double d2J0(double J, double M); // <cos^2> delta J = 0
double d2j2(double J, double M); // <cos^2> delta J = -2
double d1J1(double J, double M); // <cos> delta J = +1
double d1j1(double J, double M); // <cos> delta J = -1

double Erot(double J); // rotational energy
double E1w(double t); // electric field
double E2w(double t); // electric field

complex<double> RK(double t, int J, int M, complex<double> Cj3, complex<double> Cj2, complex<double> Cj1, complex<double> C0, complex<double> CJ1, complex<double> CJ2, complex<double> CJ3); // 4th order Runge-Kutta


int main()
{
    /**parameters and so on**/
    int Jmax = 75; //maximum value of J  < NUM
    int Jcalc;// Jcalc is the maximum value of J in calculation
    double thr = 0.0001;   // threshold of rotationl distribution
    int outnum = 0; // use in Jcaalc calculation
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
    complex<double> c[NUM], k[4][NUM]; //c[] is coefficient of wavefnction.If you change Jmax,you should change size of c[](Jmax/2+1).   k[RK][level]
    double cfin[NUM];
    complex<double> cj3 = (0.0, 0.0), cj2 = (0.0, 0.0), cj1 = (0.0, 0.0), cJ0 = (0.0, 0.0), cJ1 = (0.0, 0.0), cJ2 = (0.0, 0.0), cJ3 = (0.0, 0.0);
    complex<double> kj3 = (0.0, 0.0), kj2 = (0.0, 0.0), kj1 = (0.0, 0.0), kJ0 = (0.0, 0.0), kJ1 = (0.0, 0.0), kJ2 = (0.0, 0.0), kJ3 = (0.0, 0.0);


    /*******************************************************************************/
    double T = 0.8; // temperature in K unit
                    /*******************************************************************************/

                    /*************************************************/
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

    printf("*******************************************************************\n");
    printf(" %s", ctime(&timer));
    printf("*******************************************************************\n");
    printf(" INTENSITY (1st): %f TW \n", intensity0*pow(10.0, -12.0));
    printf(" INTENSITY (2nd): omega1 = %f TW , omega2 = %f TW\n", intensity1*pow(10.0, -12.0), rat*intensity1*pow(10.0, -12.0));
    printf(" PULSE DURATION : FWHM = %f fs\n", FWHM*pow(10.0, 15.0));
    printf(" RELATIVE PHASE : phi = %f*pi\n", phase);
    //    printf(" DELAY          : %f ps (%f*Trot)\n", delay*pow(10.0,12.0) , n );


    FILE *fp1, *fp2, *fp3, *fp4;

    fp1 = fopen("cos2.txt", "w");  // make a text file to write the final results
    fp2 = fopen("cos.txt", "w");
    fp3 = fopen("output.txt", "w");
    fp4 = fopen("pulseenvelope.txt", "w");


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


    /**calculation of partition function**/
    double b = 1 / (kB*T); // Boltzmann factor
    double Z = 0.0, z;  // Z : partition function
    double P; // statistic weight

    if (T == 0.0)
        Z = 1.0;


    else // finite temperature
    {
        for (int J = 0; J <= Jmax; J++)
        {
            z = (2.0 * double(J) + 1.0) * exp(-b * Erot(double(J)));
            Z = z + Z;
        }
    } //end of the Z calculation at the finite temperature
      // end of the calculation of partition function Z



      /**calculation of Jcalc and initial population**/

    FILE *intpop, *finpop;

    intpop = fopen("intpop.txt", "w");
    finpop = fopen("finpop.txt", "w");

    if (T == 0.0)
    {
        Jcalc = 0; // J=0 only
        for (int J = 0; J <= Jmax; J++) // J=0 only
        {
            if (J == 0)
                fprintf(intpop, "%d\t%f\n", J, 1.0);
            else
                fprintf(intpop, "%d\t%f\n", J, 0.0);
        }
    }

    else // finite temperature
    {
        for (int J = 0; J <= Jmax; J++)
        {
            z = (2.0 * double(J) + 1.0) * exp(-b * Erot(double(J)));

            if (z / Z <= thr && outnum == 0) // if the distribution can be neglected( < 100Pthr%), then set J as Jcalc
            {
                fprintf(intpop, "%d\t%f\n", J, z / Z);
                Jcalc = J - 1;
                outnum = 1;
            }

            else
            {
                fprintf(intpop, "%d\t%f\n", J, z / Z);
            }
        }
    }// end of the Jcalc calculation at the finite temperature
     // end of the calculation of Jcalc and initial population

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

            /**initial condition**/
            // thermal weight
            if (T == 0.0)
                P = 1.0;

            else // finite temperature
            {
                if (M == 0)
                    P = exp(-b * Erot(double(Jint))) / Z;
                else
                    P = 2.0 * exp(-b * Erot(double(Jint))) / Z;
            }

            N = 0; // set the number of steps 0
            step = 1;

            //initialization of the population
            for (int j = 0; j < NUM; j++)
            {
                //initial population
                if (j == Jint)
                    c[j] = complex<double>(1.0, 0.0);
                else
                    c[j] = complex<double>(0.0, 0.0);
            }
            //time evolution

            for (t = tmin; t <= tmax; t = t + dt)
            {
                //initialization
                complex<double> align = (0.0, 0.0), asum = (0.0, 0.0); // align: alignment parameter
                complex<double> orient = (0.0, 0.0), osum = (0.0, 0.0); // orient: orientaion parameter
                NORM = 0.0, Imcos2 = 0.0, Imcos = 0.0;
                dt = dtref; // small step

                if (Ethr < (E1w(t) + E2w(t))) // Runge-Kutta calculation in the region that the effect of laser pulse is important
                {

                    // calculation of k[0][]
                    for (int j = 0; j <= Jmax; j++)
                    {
                        cj3 = c[j - 3];
                        cj2 = c[j - 2];
                        cj1 = c[j - 1];
                        cJ0 = c[j];
                        cJ1 = c[j + 1];
                        cJ2 = c[j + 2];
                        cJ3 = c[j + 3];

                        // calculation of k[0][]
                        if (j < abs(M)) // impossible
                            k[0][j] = (0.0, 0.0);

                        else if (2 < j)
                        {
                            if (j == Jmax)
                            {
                                cJ1 = (0.0, 0.0);
                                cJ2 = (0.0, 0.0);
                                cJ3 = (0.0, 0.0);
                            }

                            else if (j == (Jmax - 1))
                            {
                                cJ2 = (0.0, 0.0);
                                cJ3 = (0.0, 0.0);
                            }


                            else if (j == (Jmax - 2))
                                cJ3 = (0.0, 0.0);

                            else;


                            if (j == abs(M))
                            {
                                cj1 = (0.0, 0.0);
                                cj2 = (0.0, 0.0);
                                cj3 = (0.0, 0.0);
                            }

                            else if (j == (abs(M) + 1))
                            {
                                cj2 = (0.0, 0.0);
                                cj3 = (0.0, 0.0);
                            }

                            else if (j == (abs(M) + 2))
                                cj3 = (0.0, 0.0);

                            else;

                            k[0][j] = dt * RK(t, j, M, cj3, cj2, cj1, cJ0, cJ1, cJ2, cJ3);

                        } // end of the k[0][j] calculation for 2<j


                        else if (j <= 2)
                        {
                            if (j == abs(M))
                            {
                                cj1 = (0.0, 0.0);
                                cj2 = (0.0, 0.0);
                            }

                            else if (j == (abs(M) + 1))
                                cj2 = (0.0, 0.0);

                            if (j == 2)
                                k[0][j] = dt * RK(t, j, M, 0.0, cj2, cj1, cJ0, cJ1, cJ2, cJ3);

                            else if (j == 1)
                                k[0][j] = dt * RK(t, j, M, 0.0, 0.0, cj1, cJ0, cJ1, cJ2, cJ3);

                            else if (j == 0)
                                k[0][j] = dt * RK(t, j, M, 0.0, 0.0, 0.0, cJ0, cJ1, cJ2, cJ3);

                        }    // end of the k[0][j] calculation for j <= 2


                    }// end of the k[0][j] calculation


                     // calculation of k[1][]
                    for (int j = 0; j <= Jmax; j++)
                    {
                        cj3 = c[j - 3];
                        cj2 = c[j - 2];
                        cj1 = c[j - 1];
                        cJ0 = c[j];
                        cJ1 = c[j + 1];
                        cJ2 = c[j + 2];
                        cJ3 = c[j + 3];
                        kj3 = k[0][j - 3];
                        kj2 = k[0][j - 2];
                        kj1 = k[0][j - 1];
                        kJ0 = k[0][j];
                        kJ1 = k[0][j + 1];
                        kJ2 = k[0][j + 2];
                        kJ3 = k[0][j + 3];

                        if (j < abs(M)) // impossible
                            k[1][j] = (0.0, 0.0);

                        else if (2 < j)
                        {
                            if (j == Jmax)
                            {
                                cJ1 = (0.0, 0.0);
                                cJ2 = (0.0, 0.0);
                                cJ3 = (0.0, 0.0);
                                kJ1 = (0.0, 0.0);
                                kJ2 = (0.0, 0.0);
                                kJ3 = (0.0, 0.0);
                            }

                            else if (j == (Jmax - 1))
                            {
                                cJ2 = (0.0, 0.0);
                                cJ3 = (0.0, 0.0);
                                kJ2 = (0.0, 0.0);
                                kJ3 = (0.0, 0.0);
                            }

                            else if (j == (Jmax - 2))
                            {
                                cJ3 = (0.0, 0.0);
                                kJ3 = (0.0, 0.0);
                            }

                            else;


                            if (j == abs(M))
                            {
                                cj1 = (0.0, 0.0);
                                cj2 = (0.0, 0.0);
                                cj3 = (0.0, 0.0);
                                kj1 = (0.0, 0.0);
                                kj2 = (0.0, 0.0);
                                kj3 = (0.0, 0.0);
                            }


                            else if (j == (abs(M) + 1))
                            {
                                cj2 = (0.0, 0.0);
                                cj3 = (0.0, 0.0);
                                kj2 = (0.0, 0.0);
                                kj3 = (0.0, 0.0);
                            }


                            else if (j == (abs(M) + 2))
                            {
                                cj3 = (0.0, 0.0);
                                kj3 = (0.0, 0.0);
                            }

                            else;

                            k[1][j] = dt * RK(t + 0.5 * dt, j, M, cj3 + kj3 * 0.5, cj2 + kj2 * 0.5, cj1 + kj1 * 0.5, cJ0 + kJ0 * 0.5, cJ1 + kJ1 * 0.5, cJ2 + kJ2 * 0.5, cJ3 + kJ3 * 0.5);

                        } // end of the k[1][j] calculation for 2<j


                        else if (j <= 2)
                        {

                            if (j == abs(M))
                            {
                                cj1 = (0.0, 0.0);
                                cj2 = (0.0, 0.0);
                                kj1 = (0.0, 0.0);
                                kj2 = (0.0, 0.0);
                            }

                            else if (j == (abs(M) + 1))
                            {
                                cj2 = (0.0, 0.0);
                                kj2 = (0.0, 0.0);
                            }

                            if (j == 2)
                                k[1][j] = dt * RK(t + 0.5 * dt, j, M, 0.0, cj2 + kj2 * 0.5, cj1 + kj1 * 0.5, cJ0 + kJ0 * 0.5, cJ1 + kJ1 * 0.5, cJ2 + kJ2 * 0.5, cJ3 + kJ3 * 0.5);

                            else if (j == 1)
                                k[1][j] = dt * RK(t + 0.5 * dt, j, M, 0.0, 0.0, cj1 + kj1 * 0.5, cJ0 + kJ0 * 0.5, cJ1 + kJ1 * 0.5, cJ2 + kJ2 * 0.5, cJ3 + kJ3 * 0.5);

                            else if (j == 0)
                                k[1][j] = dt * RK(t + 0.5 * dt, j, M, 0.0, 0.0, 0.0, cJ0 + kJ0 * 0.5, cJ1 + kJ1 * 0.5, cJ2 + kJ2 * 0.5, cJ3 + kJ3 * 0.5);
                        }    // end of the k[1][j] calculation for j <= 2
                    }// end of the k[1][j] calculation

                     // calculation of k[2][]
                    for (int j = 0; j <= Jmax; j++)
                    {
                        cj3 = c[j - 3];
                        cj2 = c[j - 2];
                        cj1 = c[j - 1];
                        cJ0 = c[j];
                        cJ1 = c[j + 1];
                        cJ2 = c[j + 2];
                        cJ3 = c[j + 3];

                        kj3 = k[1][j - 3];
                        kj2 = k[1][j - 2];
                        kj1 = k[1][j - 1];
                        kJ0 = k[1][j];
                        kJ1 = k[1][j + 1];
                        kJ2 = k[1][j + 2];
                        kJ3 = k[1][j + 3];

                        if (j < abs(M)) // impossible
                            k[2][j] = (0.0, 0.0);

                        else if (2 < j)
                        {
                            if (j == Jmax)
                            {
                                cJ1 = (0.0, 0.0);
                                cJ2 = (0.0, 0.0);
                                cJ3 = (0.0, 0.0);
                                kJ1 = (0.0, 0.0);
                                kJ2 = (0.0, 0.0);
                                kJ3 = (0.0, 0.0);
                            }



                            else if (j == (Jmax - 1))
                            {
                                cJ2 = (0.0, 0.0);
                                cJ3 = (0.0, 0.0);
                                kJ2 = (0.0, 0.0);
                                kJ3 = (0.0, 0.0);
                            }

                            else if (j == (Jmax - 2))
                            {
                                cJ3 = (0.0, 0.0);
                                kJ3 = (0.0, 0.0);
                            }

                            else;

                            if (j == abs(M))
                            {
                                cj1 = (0.0, 0.0);
                                cj2 = (0.0, 0.0);
                                cj3 = (0.0, 0.0);
                                kj1 = (0.0, 0.0);
                                kj2 = (0.0, 0.0);
                                kj3 = (0.0, 0.0);
                            }

                            else if (j == (abs(M) + 1))
                            {
                                cj2 = (0.0, 0.0);
                                cj3 = (0.0, 0.0);
                                kj2 = (0.0, 0.0);
                                kj3 = (0.0, 0.0);
                            }


                            else if (j == (abs(M) + 2))
                            {
                                cj3 = (0.0, 0.0);
                                kj3 = (0.0, 0.0);
                            }

                            else;

                            k[2][j] = dt * RK(t + 0.5 * dt, j, M, cj3 + kj3 * 0.5, cj2 + kj2 * 0.5, cj1 + kj1 * 0.5, cJ0 + kJ0 * 0.5, cJ1 + kJ1 * 0.5, cJ2 + kJ2 * 0.5, cJ3 + kJ3 * 0.5);

                        } // end of the k[2][j] calculation for 2<j


                        else if (j <= 2)
                        {
                            if (j == abs(M))
                            {
                                cj1 = (0.0, 0.0);
                                cj2 = (0.0, 0.0);
                                kj1 = (0.0, 0.0);
                                kj2 = (0.0, 0.0);
                            }


                            else if (j == (abs(M) + 1))
                            {
                                cj2 = (0.0, 0.0);
                                kj2 = (0.0, 0.0);
                            }


                            if (j == 2)
                                k[2][j] = dt * RK(t + 0.5 * dt, j, M, 0.0, cj2 + kj2 * 0.5, cj1 + kj1 * 0.5, cJ0 + kJ0 * 0.5, cJ1 + kJ1 * 0.5, cJ2 + kJ2 * 0.5, cJ3 + kJ3 * 0.5);

                            else if (j == 1)
                                k[2][j] = dt * RK(t + 0.5 * dt, j, M, 0.0, 0.0, cj1 + kj1 * 0.5, cJ0 + kJ0 * 0.5, cJ1 + kJ1 * 0.5, cJ2 + kJ2 * 0.5, cJ3 + kJ3 * 0.5);

                            else if (j == 0)
                                k[2][j] = dt * RK(t + 0.5 * dt, j, M, 0.0, 0.0, 0.0, cJ0 + kJ0 * 0.5, cJ1 + kJ1 * 0.5, cJ2 + kJ2 * 0.5, cJ3 + kJ3 * 0.5);

                        }    // end of the k[2][j] calculation for j <= 2

                    }// end of the k[2][j] calculation


                     // calculation of k[3][]
                    for (int j = 0; j <= Jmax; j++)
                    {
                        cj3 = c[j - 3];
                        cj2 = c[j - 2];
                        cj1 = c[j - 1];
                        cJ0 = c[j];
                        cJ1 = c[j + 1];
                        cJ2 = c[j + 2];
                        cJ3 = c[j + 3];

                        kj3 = k[2][j - 3];
                        kj2 = k[2][j - 2];
                        kj1 = k[2][j - 1];
                        kJ0 = k[2][j];
                        kJ1 = k[2][j + 1];
                        kJ2 = k[2][j + 2];
                        kJ3 = k[2][j + 3];


                        if (j < abs(M)) // impossible
                            k[2][j] = (0.0, 0.0);

                        else if (2 < j)
                        {
                            if (j == Jmax)
                            {
                                cJ1 = (0.0, 0.0);
                                cJ2 = (0.0, 0.0);
                                cJ3 = (0.0, 0.0);
                                kJ1 = (0.0, 0.0);
                                kJ2 = (0.0, 0.0);
                                kJ3 = (0.0, 0.0);
                            }

                            else if (j == (Jmax - 1))
                            {
                                cJ2 = (0.0, 0.0);
                                cJ3 = (0.0, 0.0);
                                kJ2 = (0.0, 0.0);
                                kJ3 = (0.0, 0.0);
                            }

                            else if (j == (Jmax - 2))
                            {
                                cJ3 = (0.0, 0.0);
                                kJ3 = (0.0, 0.0);
                            }

                            else;

                            if (j == abs(M))
                            {
                                cj1 = (0.0, 0.0);
                                cj2 = (0.0, 0.0);
                                cj3 = (0.0, 0.0);
                                kj1 = (0.0, 0.0);
                                kj2 = (0.0, 0.0);
                                kj3 = (0.0, 0.0);
                            }


                            else if (j == (abs(M) + 1))
                            {
                                cj2 = (0.0, 0.0);
                                cj3 = (0.0, 0.0);
                                kj2 = (0.0, 0.0);
                                kj3 = (0.0, 0.0);
                            }

                            else if (j == (abs(M) + 2))
                            {
                                cj3 = (0.0, 0.0);
                                kj3 = (0.0, 0.0);
                            }

                            else;

                            k[3][j] = dt * RK(t + dt, j, M, cj3 + kj3, cj2 + kj2, cj1 + kj1, cJ0 + kJ0, cJ1 + kJ1, cJ2 + kJ2, cJ3 + kJ3);

                        } // end of the k[3][j] calculation for 2<j


                        else if (j <= 2)
                        {
                            if (j == abs(M))
                            {
                                cj1 = (0.0, 0.0);
                                cj2 = (0.0, 0.0);
                                kj1 = (0.0, 0.0);
                                kj2 = (0.0, 0.0);
                            }


                            else if (j == (abs(M) + 1))
                            {
                                cj2 = (0.0, 0.0);
                                kj2 = (0.0, 0.0);
                            }



                            if (j == 2)
                                k[3][j] = dt * RK(t + dt, j, M, 0.0, cj2 + kj2, cj1 + kj1, cJ0 + kJ0, cJ1 + kJ1, cJ2 + kJ2, cJ3 + kJ3);

                            else if (j == 1)
                                k[3][j] = dt * RK(t + dt, j, M, 0.0, 0.0, cj1 + kj1, cJ0 + kJ0, cJ1 + kJ1, cJ2 + kJ2, cJ3 + kJ3);

                            else if (j == 0)
                                k[3][j] = dt * RK(t + dt, j, M, 0.0, 0.0, 0.0, cJ0 + kJ0, cJ1 + kJ1, cJ2 + kJ2, cJ3 + kJ3);

                        }    // end of the k[3][j] calculation for j <= 2

                    }// end of the k[3][j] calculation


                    for (int j = 0; j <= Jmax; j++)
                    {
                        c[j] = c[j] + (k[0][j] + 2.0*k[1][j] + 2.0*k[2][j] + k[3][j]) / 6.0;
                    }

                } // end of the RK calculation


                else // long duration step without RK calculation in the region that the laser pulse does not exist
                    dt = dtlong;



                //<cos^2theta> calculation
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

                //<costheta> calculation
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
                    printf("%f\t%f\t%f\t%f\t%f\n", (t + dt)*pow(10.0, 15.0), real(orient), abs(Imcos), abs(Imcos2), NORM);  // orientation, alignment check
                                                                                                                            //                    printf("%f\t%f\t%f\t%f\t%f\n",(t+dt)*pow(10.0,15.0), real(align) , real(orient) , abs(Imcos2) , NORM );  // orientation, alignment check
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

    return 0;

}


double Erot(double J)  // rotational energy in Joule unit
{
    double K = J * (J + 1.0);
    //    double ret=B*K*h*c*100.0;   // rigid rotator
    double ret = (B - D * K)*K*h*vc*100.0;   //centrifugal distortion is included.

    return ret;
}

complex<double> RK(double t, int J, int M, complex<double> Cj3, complex<double> Cj2, complex<double> Cj1, complex<double> C0, complex<double> CJ1, complex<double> CJ2, complex<double> CJ3)
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

    complex<double> ret = (I / hbar) * (term0 * C0 + term1 * CJ3 * exp(-I * wJ3*t) + term2 * Cj3 * exp(-I * wj3*t) + term3 * CJ2 * exp(-I * wJ2*t) + term4 * Cj2 * exp(-I * wj2*t) + term5 * CJ1 * exp(-I * wJ1*t) + term6 * Cj1 * exp(-I * wj1*t));
    return  ret;
}
