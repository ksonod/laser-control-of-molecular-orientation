#ifndef rotational_energy_h
#define rotational_energy_h


double Erot(double J)  // rotational energy in Joule unit
{
    double K = J * (J + 1.0);
    //    double ret=B*K*h*c*100.0;   // rigid rotator
    double ret = (B - D * K)*K*h*vc*100.0;   //centrifugal distortion is included.

    return ret;
}


double partition_function(double T)
{
    double b = 1 / (kB*T); // Boltzmann factor
    double Z = 0.0, z;  // Z : partition function
    
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

    return Z;
}


double calculate_Mlevel_population(double T, int J, int M)
{
    double P;
    double b = 1 / (kB*T); // Boltzmann factor
    double Z = partition_function(T);
    
    if (T == 0.0)
        P = 1.0;

    else // finite temperature
    {
        if (M == 0)
            P = exp(-b * Erot(double(J))) / Z;
        else
            P = 2.0 * exp(-b * Erot(double(J))) / Z;
    }
    
    return P;
}


int calculate_initial_population(double T)
{
    int Jcalc;   
    double z=0.0;
    double Z = partition_function(T);
    int outnum = 0; // use in Jcaalc calculation
    double b = 1 / (kB*T); // Boltzmann factor
    
    /*TODO: Move this parameter to somewhere*/
    double thr = 0.0001;   // threshold of rotationl distribution
 
    FILE *intpop;
    intpop = fopen("intpop.txt", "w");
    
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

    fclose(intpop);
    
    return Jcalc;
}


#endif /* rotational_energy_h */
