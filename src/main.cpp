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


int main(){
    int Jcalc;  // Jcalc is the maximum value of J in calculation
    int time_idx;  // Number of steps
    int data_sampling_counts;  // Steps in time evolution
    double rot_level_weight;  // Statistic weight
    double cos[num_time_series_data] = {0};  // time series of <costheta>
    double cos2[num_time_series_data] = {0}; // time series of <cos2theta>
    double norm_wp = 0.0;  // Norm of wavepacket coefficients.
    double calctime;  // calculation time
    double cfin[num_rot_levels] = {0};
    double t;  // Time
    double dt;  // Time step
    std::complex<double> c[num_rot_levels], k[4][num_rot_levels];  //c[] is coefficient of wavefunction.
    std::complex<double> exp_cos2; // Expectation value of cos2 (alignment parameter)
    std::complex<double> exp_cos; // Expectation value of cos (orientation parameter)
    clock_t start, end;  // calculation time

    time_t timer;
    time(&timer);

    check_params();
    Jcalc = calculate_initial_population(T);

    std::cout << "*******************************************************************\n";
    std::cout << ctime(&timer);
    std::cout << "*******************************************************************\n";
    std::cout << " INTENSITY (1st): " << intensity0*pow(10.0, -12.0) << " TW \n";
    std::cout << " INTENSITY (2nd): omega1 = " << intensity1*pow(10.0, -12.0) << " TW , omega2 = " << rat*intensity1*pow(10.0, -12.0) << " TW\n";
    std::cout << " PULSE DURATION : FWHM = " << pulse_fwhm*pow(10.0, 15.0) << " fs\n";
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

    // beginning of the calculation
    start = clock();

    // J-level calculation
    for (int Jint = 0; Jint <= Jcalc; Jint++){
        // M-level calculation
        for (int M = 0; M <= Jint; M++){
            rot_level_weight = calculate_Mlevel_population(T, Jint, M);
            
            time_idx = 0;  // set the number of steps 0
            data_sampling_counts = 0;

            //initialization of the population
            for (int j = 0; j < num_rot_levels; j++){
                //initial population
                if (j == Jint)
                    c[j] = std::complex<double>(1.0, 0.0);
                else
                    c[j] = std::complex<double>(0.0, 0.0);
            }
            
            //time evolution
            for (t = tmin; t <= tmax; t += dt){
                //initialization
                norm_wp = 0.0;
                dt = dt_small;  // small step

                if (electric_field_thr < (E1w(t) + E2w(t))){  // Runge-Kutta calculation in the region where the effect of laser pulse is important
                    // calculation of k[0][]
                    for (int j = 0; j <= Jmax; j++){
                        k[0][j] = dt * calc_schrodinger_equation(t, j, M, coef(c, j-3, M), coef(c, j-2, M), coef(c, j-1, M), coef(c, j, M), coef(c, j+1, M), coef(c, j+2, M), coef(c, j+3, M));
                    }

                    // calculation of k[1][]
                    for (int j = 0; j <= Jmax; j++){
                        k[1][j] = dt * calc_schrodinger_equation(t + 0.5 * dt, j, M, coef(c, j-3, M) + coef(k[0], j-3, M) * 0.5, coef(c, j-2, M) + coef(k[0], j-2, M) * 0.5, coef(c, j-1, M) + coef(k[0], j-1, M) * 0.5, coef(c, j, M) + coef(k[0], j, M) * 0.5, coef(c, j+1, M) + coef(k[0], j+1, M) * 0.5, coef(c, j+2, M) + coef(k[0], j+2, M) * 0.5, coef(c, j+3, M) + coef(k[0], j+3, M) * 0.5);
                    }

                    // calculation of k[2][]
                    for (int j = 0; j <= Jmax; j++){
                        k[2][j] = dt * calc_schrodinger_equation(t + 0.5 * dt, j, M, coef(c, j-3, M) + coef(k[1], j-3, M) * 0.5, coef(c, j-2, M) + coef(k[1], j-2, M) * 0.5, coef(c, j-1, M) + coef(k[1], j-1, M) * 0.5, coef(c, j, M) + coef(k[1], j, M) * 0.5, coef(c, j+1, M) + coef(k[1], j+1, M) * 0.5, coef(c, j+2, M) + coef(k[1], j+2, M) * 0.5, coef(c, j+3, M) + coef(k[1], j+3, M) * 0.5);
                    }

                    // calculation of k[3][]
                    for (int j = 0; j <= Jmax; j++){
                        k[3][j] = dt * calc_schrodinger_equation(t + dt, j, M, coef(c, j-3, M) + coef(k[2], j-3, M), coef(c, j-2, M) + coef(k[2], j-2, M), coef(c, j-1, M) + coef(k[2], j-1, M), coef(c, j, M) + coef(k[2], j, M), coef(c, j+1, M) + coef(k[2], j+1, M), coef(c, j+2, M) + coef(k[2], j+2, M), coef(c, j+3, M) + coef(k[2], j+3, M));
                    }

                    for (int j = 0; j <= Jmax; j++){
                        c[j] += (k[0][j] + 2.0*k[1][j] + 2.0*k[2][j] + k[3][j]) / 6.0;
                    }
                    
                    data_sampling_counts +=1 ;  // Data will be sampled with the interval of data_sampling_step.
                }  // end of the RK calculation
                else{  // Large step without RK calculation in the region where the laser pulse is extremely weak or does not exist.
                    dt = dt_large;
                    data_sampling_counts = 0;  // Data will be always sampled because of the large temporal step.
                }

                exp_cos2 = calculate_cos2_expectation_value(c, M, t+dt);  // Expectation value of cos2
                exp_cos = calculate_cos_expectation_value(c, M, t+dt);  // Expectation value of cos

                if (T == 0.0){
                    for (int j = 0; j <= Jmax; j++){
                        norm_wp += norm(c[j]);  // Norm calculation
                    }
                    std::cout << (t + dt)*pow(10.0, 15.0) << "\t" << exp_cos.real() << "\t" << exp_cos.imag() << "\t" << exp_cos2.imag() << "\t" << norm_wp << "\n";
                }
                      
                // Storing and saving expectation values of cos2 and cos.
                if (data_sampling_counts%data_sampling_step == 0){
                    cos2[time_idx] += rot_level_weight * exp_cos2.real();  // Time series of <cos^2theta>
                    cos[time_idx] += rot_level_weight * exp_cos.real();  // Time series of <costheta>

                    if (Jint == Jcalc && M == Jint){  // Write results in the text file at the end of J & M loop.
                        file_cos2 << (t + dt) / t_rot_period << "\t" << cos2[time_idx] << "\n";
                        file_cos << (t + dt) / t_rot_period << "\t" << cos[time_idx] << "\n";
                    }
                    time_idx += 1; //next step data
                }
            }  // end of the time evolution

             // Calculation of norm and final population
            for (int j = 0; j <= Jmax; j++){
                if (0.0 < T){
                    norm_wp += norm(c[j]);  // Norm calculation
                }
                
                cfin[j] += rot_level_weight * norm(c[j]); // Final population
                if (Jint == Jcalc && M == Jint){  // Write the results in txt file at the end
                    file_finpop << j << "\t" << cfin[j] << "\n";
                }
            }  // end of the final population calculation

            std::cout << "J" << Jint << "M" << M << "\tN=" << norm_wp << "\tIm(cos)=" << exp_cos.imag() << "\tIm(cos2)=" << exp_cos2.imag() << "\n";
            file_output << " J" << Jint << "M" << M << "\tN=" << norm_wp << "\tIm(cos)=" << exp_cos.imag() << "\tIm(cos2)=" << exp_cos2.imag() << "\n";
        }  // end of the calculation of all M for specific J
    }  // end of the all Jint calculation

    end = clock(); // calculation time
    calctime = (double(end - start)) / CLOCKS_PER_SEC; // calculation time (seconds)

    auto [cosmax, cosmin] = return_max_min_values(cos);
    auto [cos2max, cos2min] = return_max_min_values(cos2);

    // output the maximum value of <costheta>
    std::cout << "-----------------------\n";
    std::cout << "calculation time = " << calctime << " sec.\n";
    std::cout << "<cos>max  = " << cosmax << "\n";
    std::cout << "<cos>min  = " << cosmin << "\n";
    std::cout << "<cos2>max = " << cos2max << "\n";

    // write the conditions and results in txt file
    file_output << "*******************************************************************\n";
    file_output << "-- LASER PULSE --\n";
    file_output << " -first pulse\n";
    file_output << " INTENSITY      : " << intensity0*pow(10.0, -12.0) << " TW\n";
    file_output << " PULSE DURATION : FWHM = " << pulse_fwhm*pow(10.0, 15.0) << " fs\n";
    file_output << " -second pulse\n";
    file_output << " INTENSITY      : omega1 = " << intensity1*pow(10.0, -12.0) << " TW , omega2 = " << rat*intensity1*pow(10.0, -12.0) << " TW\n";
    file_output << " PULSE DURATION : FWHM = " << pulse_fwhm*pow(10.0, 15.0) << " fs\n";
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
