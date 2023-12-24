#ifndef matrix_element_h
#define matrix_element_h


double d1j1(double J, double M) // <JM|cos|J-1M>
{
    if (J < abs(M) || J == 0){
        return 0.0;
    }
    else{
        double a = (J + M)*(J - M) / (2.0 * J + 1.0) / (2.0 * J - 1.0);
        double ret = sqrt(a);
        return ret;
    }
}



double d1J1(double J, double M) // <JM|cos|J+1M>
{
    double a = (J + M + 1.0)*(J - M + 1.0) / (2.0*J + 3.0) / (2.0*J + 1.0);
    double ret = sqrt(a);

    if ((J + 1.0) < abs(M))
    {
        return 0.0;
    }

    else
    {
        return ret;
    }
}



double d2j2(double J, double M) // <JM|cos^2|J-2M>
{
    if ((abs(M) <= J) && (2.0 <= J))
    {
        double a = (J * J - M * M) * ((J - 1.0) * (J - 1.0) - M * M) / (2.0 * J + 1.0) / (2.0 * J - 3.0);
        double ret = 1.0 / (2.0 * J - 1.0) * sqrt(a);
        return ret;
    }

    else
        return 0.0;
}


double d2J0(double J, double M)  // <JM|cos^2|JM>
{
    if (abs(M) <= J)
    {
        double a = (J * (J + 1.0) - 3.0 * M * M) / (2.0 * J + 3.0) / (2.0 * J - 1.0);
        double ret = 2.0 / 3.0 * a + 1.0 / 3.0;
        return ret;
    }

    else
        return 0.0;
}


double d2J2(double J, double M)  // <JM|cos^2|J+2M>
{
    if (abs(M) <= J)
    {
        double a = ((J + 2.0) * (J + 2.0) - M * M) * ((J + 1.0) * (J + 1.0) - M * M) / (2.0 * J + 5.0) / (2.0 * J + 1.0);
        double ret = 1.0 / (2.0 * J + 3.0) * sqrt(a);

        return ret;
    }

    else
        return 0.0;
}


double d3j3(double J, double M) // <JM|cos3|J-3M>
{

    if ((abs(M) <= J) && (3.0 <= J))
    {
        double a = (J + M)*(J + M - 1.0)*(J + M - 2.0)*(J - M)*(J - M - 1.0)*(J - M - 2.0) / (2.0*J + 1.0) / (2.0*J - 5.0);
        double ret = sqrt(a) / (2.0 * J - 1.0) / (2.0 * J - 3.0);
        return ret;
    }

    else
        return 0.0;
}


double d3j1(double J, double M) // <JM|cos3|J-1M>
{
    if ((abs(M) <= J) && (1.0 <= J))
    {
        double a = (J + M)*(J - M) / (2.0 * J + 1.0) / (2.0 * J - 1.0);
        double b = 3.0 * (J * J - M * M - 2.0) / (2.0 * J + 3.0) / (2.0 * J - 3.0);
        double ret = b * sqrt(a);
        return ret;
    }

    else
        return 0.0;
}


double d3J1(double J, double M) // <JM|cos3|J+1M>
{
    if (abs(M) <= J)
    {
        double a = (J + M + 1.0)*(J - M + 1.0) / (2.0 * J + 3.0) / (2.0 * J + 1.0);
        double b = 3.0 * (J * (J + 2.0) - M * M - 1.0) / (2.0 * J + 5.0) / (2.0 * J - 1.0);
        double ret = b * sqrt(a);
        return ret;
    }

    else
        return 0.0;
}


double d3J3(double J, double M) // <JM|cos3|J+3M>
{
    if (abs(M) <= J)
    {
        double a = (J - M + 1.0)*(J - M + 2.0)*(J - M + 3.0)*(J + M + 1.0)*(J + M + 2.0)*(J + M + 3.0) / (2.0 * J + 7.0) / (2.0 * J + 1.0);
        double ret = sqrt(a) / (2.0 * J + 5.0) / (2.0 * J + 3.0);

        return ret;
    }

    else
        return 0.0;
}


#endif /* matrix_element_h */
