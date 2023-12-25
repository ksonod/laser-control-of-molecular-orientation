#ifndef rotational_energy_h
#define rotational_energy_h


// Rotational energy in Joule unit
double rot_energy(int j){
    double K = double(j) * (double(j) + 1.0);
    double ret = (B_rot - D_centr * K) * K * h_planck * speed_of_light * 100.0;   // Centrifugal distortion is included.
    return ret;
}


double partition_function(double T){
    double Z = 0.0, z;  // Z : partition function
    
    if (T == 0.0){
        Z = 1.0;
    }
    else{
        double b = 1 / (k_boltzmann*T); // Boltzmann factor

        for (int j = 0; j <= Jmax; j++){
            z = (2.0 * double(j) + 1.0) * exp(-b * rot_energy(j));
            Z = z + Z;
        }
    }
    return Z;
}


double calculate_Mlevel_population(double T, int j, int m){
    double p;
    
    if (T == 0.0){
        p = 1.0;
    }
    else{ // finite temperature
        double b = 1 / (k_boltzmann*T); // Boltzmann factor
        double Z = partition_function(T);

        if (m == 0){
            p = exp(-b * rot_energy(j)) / Z;
        }
        else{
            p = 2.0 * exp(-b * rot_energy(j)) / Z;
        }
    }
    
    return p;
}


int calculate_initial_population(double T){
    int Jcalc = 0 ;
    double z=0.0;
    double Z = partition_function(T);
    int outnum = 0; // use in Jcaalc calculation
     
    FILE* intpop = std::fopen("intpop.txt", "w");
    
    if (T == 0.0){
        Jcalc = 0; // J=0 only
        for (int j = 0; j <= Jmax; j++){ // J=0 only
            if (j == 0){
                fprintf(intpop, "%d\t%f\n", j, 1.0);
            }
            else{
                fprintf(intpop, "%d\t%f\n", j, 0.0);
            }
        }
    }
    else{  // finite temperature
        double b = 1 / (k_boltzmann*T); // Boltzmann factor

        for (int j = 0; j <= Jmax; j++){
            z = (2.0 * double(j) + 1.0) * exp(-b * rot_energy(j));

            if (z / Z <= rot_population_thr && outnum == 0){  // if the distribution can be neglected( < 100Pthr%), then set J as Jcalc
                fprintf(intpop, "%d\t%f\n", j, z / Z);
                Jcalc = j - 1;
                outnum = 1;
            }
            else{
                fprintf(intpop, "%d\t%f\n", j, z / Z);
            }
        }
    }

    std::fclose(intpop);
    
    return Jcalc;
}

#endif /* rotational_energy_h */
