#ifndef expectation_values_h
#define expectation_values_h


std::complex<double> calculate_cos2_expectation_value(std::complex<double> (&c)[num_rot_levels], int m, double t){

    std::complex<double> exp_cos2 = 0.0;

    for (int j = 0; j <= Jmax; j++){
        exp_cos2 += conj(coef(c, j, m)) * coef(c, j, m) * d2J0(j, m) + conj(coef(c, j, m)) * coef(c, j - 2, m) * d2j2(j, m) * exp(-I * (rot_energy(j - 2) - rot_energy(j)) * t / hbar) + conj(coef(c, j, m)) * coef(c, j + 2, m) * d2J2(j, m) * exp(-I * (rot_energy(j + 2) - rot_energy(j)) * t / hbar);
    }
    
    return exp_cos2;
}


std::complex<double> calculate_cos_expectation_value(std::complex<double> (&c)[num_rot_levels], int m, double t){

    std::complex<double> exp_cos = 0.0;

    for (int j = 0; j <= Jmax; j++){
        exp_cos += conj(coef(c, j, m)) * coef(c, j - 1, m) * d1j1(j, m) * exp(-I * (rot_energy(j - 1) - rot_energy(j)) * t / hbar) + conj(coef(c, j, m)) * coef(c, j + 1, m) * d1J1(j, m) * exp(-I * (rot_energy(j + 1) - rot_energy(j)) * t / hbar);
    }

    return exp_cos;
}


std::tuple<double, double> return_max_min_values(double (&arr)[num_time_series_data]){
    
    double max_val = 0.0;
    double min_val = 0.0;
    
    for (int i=0; i<num_time_series_data; i++){
        if (max_val < arr[i]){
            max_val = arr[i];
        }
        
        if (min_val > arr[i]){
            min_val = arr[i];
        }
    }
    
    return {max_val, min_val};
}


#endif /* expectation_values_h */
