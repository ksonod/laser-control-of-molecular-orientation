#ifndef matrix_element_h
#define matrix_element_h


// <JM|cos|J-1M>
double d1j1(double j, double m){
    if (j <= abs(m) || j == 0.0){
        return 0.0;
    }
    else{
        return sqrt( (j + m)*(j - m) / (2.0 * j + 1.0) / (2.0 * j - 1.0) );
    }
}


// <JM|cos|J+1M>
double d1J1(double j, double m){
    if (j < abs(m)){
        return 0.0;
    }
    else{
        return sqrt( (j + m + 1.0)*(j - m + 1.0) / (2.0*j + 3.0) / (2.0*j + 1.0) );
    }
}


// <JM|cos^2|J-2M>
double d2j2(double j, double m){
    if ((j < 2.0) || abs(m) > j-2){
        return 0.0;
    }
    else{
        double a = (j * j - m * m) * ((j - 1.0) * (j - 1.0) - m * m) / (2.0 * j + 1.0) / (2.0 * j - 3.0);
        double ret = 1.0 / (2.0 * j - 1.0) * sqrt(a);
        return ret;
    }
}


// <JM|cos^2|JM>
double d2J0(double j, double m){
    if (j < abs(m)){
        return 0.0;
    }
    else{
        double a = (j * (j + 1.0) - 3.0 * m * m) / (2.0 * j + 3.0) / (2.0 * j - 1.0);
        double ret = 2.0 / 3.0 * a + 1.0 / 3.0;
        return ret;
    }
}


// <JM|cos^2|J+2M>
double d2J2(double j, double m){
    if (j < abs(m)){
        return 0.0;
    }
    else{
        double a = ((j + 2.0) * (j + 2.0) - m * m) * ((j + 1.0) * (j + 1.0) - m * m) / (2.0 * j + 5.0) / (2.0 * j + 1.0);
        double ret = 1.0 / (2.0 * j + 3.0) * sqrt(a);
        return ret;
    }
}


// <JM|cos3|J-3M>
double d3j3(double j, double m){
    if ((j < 3) || (j-3 < abs(m))){
        return 0.0;
    }
    else{
        double a = (j + m)*(j + m - 1.0)*(j + m - 2.0)*(j - m)*(j - m - 1.0)*(j - m - 2.0) / (2.0*j + 1.0) / (2.0*j - 5.0);
        double ret = sqrt(a) / (2.0 * j - 1.0) / (2.0 * j - 3.0);
        return ret;
    }
}


// <JM|cos3|J-1M>
double d3j1(double j, double m){
    if ((j < 1) || (j-1 < abs(m))){
        return 0.0;
    }
    else{
        double a = (j + m)*(j - m) / (2.0 * j + 1.0) / (2.0 * j - 1.0);
        double b = 3.0 * (j * j - m * m - 2.0) / (2.0 * j + 3.0) / (2.0 * j - 3.0);
        double ret = b * sqrt(a);
        return ret;
    }
}


// <JM|cos3|J+1M>
double d3J1(double j, double m){
    if (j < abs(m)){
        return 0.0;
    }
    else{
        double a = (j + m + 1.0)*(j - m + 1.0) / (2.0 * j + 3.0) / (2.0 * j + 1.0);
        double b = 3.0 * (j * (j + 2.0) - m * m - 1.0) / (2.0 * j + 5.0) / (2.0 * j - 1.0);
        double ret = b * sqrt(a);
        return ret;
    }
}


// <JM|cos3|J+3M>
double d3J3(double j, double m){
    if (j < abs(m)){
        return 0.0;
    }
    else{
        double a = (j - m + 1.0)*(j - m + 2.0)*(j - m + 3.0)*(j + m + 1.0)*(j + m + 2.0)*(j + m + 3.0) / (2.0 * j + 7.0) / (2.0 * j + 1.0);
        double ret = sqrt(a) / (2.0 * j + 5.0) / (2.0 * j + 3.0);
        return ret;
    }
}


#endif /* matrix_element_h */
