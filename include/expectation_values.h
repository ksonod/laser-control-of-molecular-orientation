#ifndef expectation_values_h
#define expectation_values_h


std::complex<double> calculate_cos2_expectation_value(std::complex<double> (&c)[num_rot_levels], int m, double t, double dt){

    std::complex<double> exp_cos2 = 0.0, temp_cos2 = 0.0;

    for (int j = 0; j <= Jmax; j++){
        if (j < abs(m)) // impossible
            temp_cos2 = 0.0;

        else if (j == Jmax || j == Jmax - 1){
            if (j == abs(m) || j == abs(m) + 1){
                temp_cos2 = conj(c[j])*c[j] * d2J0(j, m);
            }
            else{
                temp_cos2 = conj(c[j])*c[j] * d2J0(j, m) + conj(c[j])*c[j - 2] * d2j2(j, m)*exp(-I * (rot_energy(j - 2) - rot_energy(j)) * (t + dt) / hbar);
            }
            exp_cos2 += temp_cos2;
        }
        else if (j == 0 || j == 1){
            temp_cos2 = conj(c[j])*c[j] * d2J0(j, m) + conj(c[j])*c[j + 2] * d2J2(j, m)*exp(-I * (rot_energy(j + 2) - rot_energy(j)) * (t + dt) / hbar);
            exp_cos2 += temp_cos2;
        }
        else{
            if (j - 2 < abs(m)){
                temp_cos2 = conj(c[j])*c[j] * d2J0(j, m) + conj(c[j])*c[j + 2] * d2J2(j, m)*exp(-I * (rot_energy(j + 2) - rot_energy(j)) * (t + dt) / hbar);
            }
            else{
                temp_cos2 = conj(c[j])*c[j] * d2J0(j, m) + conj(c[j])*c[j - 2] * d2j2(j, m)*exp(-I * (rot_energy(j - 2) - rot_energy(j)) * (t + dt) / hbar) + conj(c[j])*c[j + 2] * d2J2(j, m)*exp(-I * (rot_energy(j + 2) - rot_energy(j)) * (t + dt) / hbar);
            }
            exp_cos2 += temp_cos2;
        }
    }
    
    return exp_cos2;
}


std::complex<double> calculate_cos_expectation_value(std::complex<double> (&c)[num_rot_levels], int m, double t, double dt){

    std::complex<double> exp_cos = 0.0, temp_cos = 0.0;

    for (int j = 0; j <= Jmax; j++){
        if (j < abs(m)){
            temp_cos = 0.0;
        }
        else if (j == Jmax){ // no interaction with upper state
            if (j == abs(m)){ // no interaction
                temp_cos = 0.0;
            }
            else{ // interaction with lower state(J=J-1)
                temp_cos = conj(c[j])*c[j - 1] * d1j1(j, m)*exp(-I * (rot_energy(j - 1) - rot_energy(j)) * (t + dt) / hbar);
            }
            exp_cos += temp_cos;
        }
        else if (j == 0){ // no interaction with lower state
            temp_cos = conj(c[j])*c[j + 1] * d1J1(j, m)*exp(-I * (rot_energy(j + 1) - rot_energy(j)) * (t + dt) / hbar);
            exp_cos += temp_cos;
        }
        else{
            if (j == abs(m)){ // no interaction with lower state
                temp_cos = conj(c[j])*c[j + 1] * d1J1(j, m)*exp(-I * (rot_energy(j + 1) - rot_energy(j)) * (t + dt) / hbar);
            }
            else{
                temp_cos = conj(c[j])*c[j - 1] * d1j1(j, m)*exp(-I * (rot_energy(j - 1) - rot_energy(j)) * (t + dt) / hbar) + conj(c[j])*c[j + 1] * d1J1(j, m)*exp(-I * (rot_energy(j + 1) - rot_energy(j)) * (t + dt) / hbar);
            }
            exp_cos += temp_cos;
        }
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
