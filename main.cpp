#include <random>
#include <iostream>
#include <ctime>
#include <fstream>
#include <vector>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>

#include "BlackScholes.hpp"


int main(int argc, const char * argv[]) {
    
    
    BlackScholes Market_test(0.05,0.2,1,100,100);
    
    
    
    // EDP
    
    /*
    float maxTemps = 1;
    float maxS = log(300);
    int nbr_T = 50000; // M dans l'énoncé
    int nbr_S = 1000; // N dans l'énoncé
    
    //Market_test.write_elem_finie(maxTemps,maxS, nbr_T, nbr_S,1);
    //Market_test.write_elem_finie(maxTemps,maxS, nbr_T, nbr_S,2);
    //Market_test.write_elem_finie(maxTemps,maxS, nbr_T, nbr_S,3);
    //Market_test.write_exact_BS(nbr_S, maxS);
    
    int nbr_T_0 = 10000;
    int pas_T = 5000;
    int nbr_val = 10;
    
    // Mesure de l'erreur par rapport à M
    
    Market_test.write_EDP_M(maxTemps,maxS, nbr_T_0, pas_T,nbr_S,1, nbr_val);
    Market_test.write_EDP_M(maxTemps,maxS, nbr_T_0, pas_T,nbr_S,2, nbr_val);
    Market_test.write_EDP_M(maxTemps,maxS, nbr_T_0, pas_T,nbr_S,3, nbr_val);
    
    
    // Mesure de l'erreur par rapport à N
    
    int nbr_S_0 = 500;
    int pas_S = 100;
     */
    //int nbr_val = 30;
    
    //Market_test.write_EDP_N(maxTemps,maxS, nbr_T, pas_S,nbr_S_0,1, nbr_val);
    //Market_test.write_EDP_N(maxTemps,maxS, nbr_T, pas_S,nbr_S_0,2, nbr_val);
    //Market_test.write_EDP_N(maxTemps,maxS, nbr_T, pas_S,nbr_S_0,3, nbr_val);
    
    
    // EDP pour American et Bermudan
    /*
    float maxTemps = 1;
    float maxS = log(300);
    int nbr_T = 50000; // M dans l'énoncé
    int nbr_S = 1000; // N dans l'énoncé
    Market_test.write_elem_finie(maxTemps,maxS, nbr_T, nbr_S,4);
    Market_test.write_elem_finie(maxTemps,maxS, nbr_T, nbr_S,5);
    */
    

    return 0;
}

