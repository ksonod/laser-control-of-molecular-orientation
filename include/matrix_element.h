#ifndef matrix_element_h
#define matrix_element_h


// <JM|cos|J-1M>
double d1j1(int j, int m){
    if (j <= abs(m) || j == 0){
        return 0.0;
    }
    else{
        double J = double(j);
        double M = double(m);
        return sqrt( (J + M) * (J - M) / (2.0 * J + 1.0) / (2.0 * J - 1.0) );
    }
}


// <JM|cos|J+1M>
double d1J1(int j, int m){
    if (j < abs(m)){
        return 0.0;
    }
    else{
        double J = double(j);
        double M = double(m);
        return sqrt( (J + M + 1.0) * (J - M + 1.0) / (2.0 * J + 3.0) / (2.0 * J + 1.0) );
    }
}


// <JM|cos^2|J-2M>
double d2j2(int j, int m){
    if ((j < 2) || j-2 < abs(m)){
        return 0.0;
    }
    else{
        double J = double(j);
        double M = double(m);
        double a = (J * J - M * M) * ((J - 1.0) * (J - 1.0) - M * M) / (2.0 * J + 1.0) / (2.0 * J - 3.0);
        return 1.0 / (2.0 * J - 1.0) * sqrt(a);
    }
}


// <JM|cos^2|JM>
double d2J0(int j, int m){
    if (j < abs(m)){
        return 0.0;
    }
    else{
        double J = double(j);
        double M = double(m);
        double a = (J * (J + 1.0) - 3.0 * M * M) / (2.0 * J + 3.0) / (2.0 * J - 1.0);
        return 2.0 / 3.0 * a + 1.0 / 3.0;
    }
}


// <JM|cos^2|J+2M>
double d2J2(int j, int m){
    if (j < abs(m)){
        return 0.0;
    }
    else{
        double J = double(j);
        double M = double(m);
        double a = ((J + 2.0) * (J + 2.0) - M * M) * ((J + 1.0) * (J + 1.0) - M * M) / (2.0 * J + 5.0) / (2.0 * J + 1.0);
        return 1.0 / (2.0 * J + 3.0) * sqrt(a);
    }
}


// <JM|cos3|J-3M>
double d3j3(int j, int m){
    if ((j < 3) || (j-3 < abs(m))){
        return 0.0;
    }
    else{
        double J = double(j);
        double M = double(m);
        double a = (J + M) * (J + M - 1.0) * (J + M - 2.0) * (J - M) * (J - M - 1.0) * (J - M - 2.0) / (2.0 * J + 1.0) / (2.0 * J - 5.0);
        return sqrt(a) / (2.0 * J - 1.0) / (2.0 * J - 3.0);
    }
}


// <JM|cos3|J-1M>
double d3j1(int j, int m){
    if ((j < 1) || (j-1 < abs(m))){
        return 0.0;
    }
    else{
        double J = double(j);
        double M = double(m);
        double a = (J + M) * (J - M) / (2.0 * J + 1.0) / (2.0 * J - 1.0);
        double b = 3.0 * (J * J - M * M - 2.0) / (2.0 * J + 3.0) / (2.0 * J - 3.0);
        return b * sqrt(a);
    }
}


// <JM|cos3|J+1M>
double d3J1(int j, int m){
    if (j < abs(m)){
        return 0.0;
    }
    else{
        double J = double(j);
        double M = double(m);
        double a = (J + M + 1.0) * (J - M + 1.0) / (2.0 * J + 3.0) / (2.0 * J + 1.0);
        double b = 3.0 * (J * (J + 2.0) - M * M - 1.0) / (2.0 * J + 5.0) / (2.0 * J - 1.0);
        return b * sqrt(a);
    }
}


// <JM|cos3|J+3M>
double d3J3(int j, int m){
    if (j < abs(m)){
        return 0.0;
    }
    else{
        double J = double(j);
        double M = double(m);
        double a = (J - M + 1.0) * (J - M + 2.0) * (J - M + 3.0) * (J + M + 1.0) * (J + M + 2.0) * (J + M + 3.0) / (2.0 * J + 7.0) / (2.0 * J + 1.0);
        return sqrt(a) / (2.0 * J + 5.0) / (2.0 * J + 3.0);
    }
}


#endif /* matrix_element_h */
