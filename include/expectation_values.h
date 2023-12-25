#ifndef expectation_values_h
#define expectation_values_h


std::complex<double> calculate_cos2_expectation_value(std::complex<double> (&c)[num_rot_levels], int M, double t, double dt){

    std::complex<double> align = 0.0, asum = 0.0; // align: alignment parameter

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
    
    return align;
}


std::complex<double> calculate_cos_expectation_value(std::complex<double> (&c)[num_rot_levels], int M, double t, double dt){

    std::complex<double> orient = 0.0, osum = 0.0; // orient: orientaion parameter

    for (int j = 0; j <= Jmax; j++){
        if (j < abs(M))
            osum = 0.0;

        else if (j == Jmax){ // no interaction with upper state
            if (j == abs(M)) // no interaction
                osum = 0.0;
            else // interaction with lower state(J=J-1)
                osum = conj(c[j])*c[j - 1] * d1j1(double(j), double(M))*exp(-I * (Erot(double(j) - 1.0) - Erot(double(j))) * (t + dt) / hbar);
            orient = osum + orient;
        }

        else if (j == 0){ // no interaction with lower state
            osum = conj(c[j])*c[j + 1] * d1J1(double(j), double(M))*exp(-I * (Erot(double(j) + 1.0) - Erot(double(j))) * (t + dt) / hbar);
            orient = osum + orient;
        }

        else{
            if (j == abs(M)) // no interaction with lower state
                osum = conj(c[j])*c[j + 1] * d1J1(double(j), double(M))*exp(-I * (Erot(double(j) + 1.0) - Erot(double(j))) * (t + dt) / hbar);
            else
                osum = conj(c[j])*c[j - 1] * d1j1(double(j), double(M))*exp(-I * (Erot(double(j) - 1.0) - Erot(double(j))) * (t + dt) / hbar) + conj(c[j])*c[j + 1] * d1J1(double(j), double(M))*exp(-I * (Erot(double(j) + 1.0) - Erot(double(j))) * (t + dt) / hbar);
            orient = osum + orient;
        }
    }// end of the <costheta> calculation

    return orient;
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
