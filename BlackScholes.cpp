#include "BlackScholes.hpp"
#include <math.h>
#include <tgmath.h>
#include <random>
#include <ctime>

#include <iostream>
#include <fstream>


BlackScholes::BlackScholes(float r, float s, float T, float S0, float K ) : r(r), sigma(s), T(T), S0(S0), K(K) {
};


float BlackScholes::get_C(){
    float x = S0*phi(get_d1());
    float y = K*exp(-1*r*T)*phi(get_d2());
    return x-y;
};


float BlackScholes::phi(float d){
    return 0.5*(1+erf(d/sqrt(2)));
};


float BlackScholes::get_d1(){
    return (log(S0/K)+(r+pow(sigma,2)/2)*T)/(sigma*sqrt(T));
};


float BlackScholes::get_d2(){
    return get_d1()- sigma * sqrt(T);
};


float BlackScholes::get_Z(){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution <> dis(0.0, 1.0);
    
    float V = dis(gen);
    float U = dis(gen);
    
    float Z = sqrt(-2*log(U))*cos(2*M_PI*V);
    return Z;
};


float BlackScholes::get_X(){
    
    float X = S0*exp((r-pow(sigma,2)/2)*T+sigma*sqrt(T)*get_Z())-K;
    if(X>0){
        return X * exp(-r*T);
    }
    else{
        return 0;
    }
};

void BlackScholes::change_S0(float s){
    S0 = s;
};


float* BlackScholes::MC(int n){
    float *array= new float[2];
    
    float Esp = 0;
    float Esp2 = 0;
    
    for(int i=0 ; i< n ; ++i){
        float X = get_X();
        Esp += X;
        Esp2 += pow(X,2);
    }
    array[0]=Esp/n;
    array[1]=Esp2/n-pow(Esp/n,2);
    return array;
};


float BlackScholes::eval_time_MC(int n, int m){
    
    struct timespec start , stop;
    clock_gettime(CLOCK_MONOTONIC ,&start);
    for(int i=0 ; i< m ; ++i){
        MC(n);
    }
    clock_gettime(CLOCK_MONOTONIC ,&stop);
    double elapsed = (stop.tv_sec - start.tv_sec);
    elapsed += (stop.tv_nsec - start.tv_nsec) / 1000000000.0;
    
    return elapsed/m;
};


float* BlackScholes::eval_val_MC(int n, int m){
    
    float *res_MC= new float[2];
    float *array= new float[2];
    float mean_MC = 0;
    float var_MC = 0;
    for(int i=0 ; i< m ; ++i){
        res_MC =MC(n);
        mean_MC += res_MC[0];
        var_MC += res_MC[1];
    }
    array[0]=mean_MC/m;
    array[1]=var_MC/m;
    
    return array;
};

void BlackScholes::write_MC(int m, int N, int k){
    std::ofstream outputFile;
    outputFile.open("write_MC.csv");
    for(int j=0;j<N;++j){
        outputFile << eval_val_MC((1+j)*k,m)[0]-get_C() << '\n' ;
        outputFile.close();
    }
};


float BlackScholes::get_next_S(float St, float pas){
    float S_next = 0;
    S_next = St * exp( (r-pow(sigma,2)/2) * pas + sigma * sqrt(pas) * get_Z() );
    return S_next;
};



std::vector<double> BlackScholes::element_finie(float maxTemps,float maxS, int nbr_T, int nbr_S,int val){
    // initialisation des variables
    std::vector<double> u(2*nbr_S+2, 0.0);
    u[0]=0;
    u[2*nbr_S] = exp(maxS)-K;
    float y=0;
    for(int j =1;j<(2*nbr_S);++j){
        y = (j-nbr_S)*maxS/(nbr_S);
        if( (exp(y)-K)>0 ){
            u[j] = exp(y) - K ;
        }
        else{
            u[j] = 0;
        }
    }
    // initialisation de A
    float h = 2*maxS/(2*nbr_S+1);
    float a = 0.5*pow(sigma/h,2)-(r-pow(sigma,2)*0.5)*0.5/h;    // valeur de gauche ie terme en "i-1"
    float b = -pow(sigma/h,2)-r;                                // valeur du milieu ie terme en "i"
    float c = 0.5*pow(sigma/h,2)+(r-pow(sigma,2)*0.5)*0.5/h;    // valeur de droite ie terme en "i+1"
    if(val==4 or val==5){
        u[0]=K-exp(maxS);
        u[2*nbr_S] = 0;
        float y=0;
        for(int j =1;j<(2*nbr_S);++j){
            y = (j-nbr_S)*maxS/(nbr_S);
            if( (K-exp(y))>0 ){
                u[j] = K-exp(y) ;
            }
            else{
                u[j] = 0;
            }
        }
    }
    
    
    float dist=0;
    float dist_nouv =0;
    // boucle sur le temps
    for(int j=1;j<nbr_T;++j){
        dist_nouv=0;
        std::vector<double> v = val_exact_BS_T(nbr_S, maxS,(j)*T/nbr_T);
        // calcul de la matrice des temps suivants
        if(val==1){
            u = next_euler_expl(u,a,b,c,2*nbr_S+1,maxTemps/nbr_T);
            u[2*nbr_S] = exp(y) - K*exp(-j*T*r/nbr_T);
        }
        if(val==2){
            u = next_euler_impl(u ,a,b,c,2*nbr_S+1,maxTemps/nbr_T);
            u[2*nbr_S] = exp(y) - K*exp(-j*T*r/nbr_T);
        }
        if(val==3){
            u = next_Crank_Nicholson(u,a,b,c,2*nbr_S+1,maxTemps/nbr_T,0.5);
            u[2*nbr_S] = exp(y) - K*exp(-j*T*r/nbr_T);
        }
        if(val==4){
            u = next_American_impl(u ,a,b,c,2*nbr_S+1,maxTemps/nbr_T,h,maxS);
            u[2*nbr_S] = 0; // car put
        }
        if(val==5){
            int nbr_Steps = 100;
            u = next_Bermudan_impl(u ,a,b,c,2*nbr_S+1,maxTemps/nbr_T,h,maxS,nbr_Steps);
            u[2*nbr_S] = 0; // car put
        }
        
        for(int j=0;j<(2*nbr_S+1);++j){
            dist_nouv = pow(u[j]-v[j],2)/(2*nbr_S+1)+dist_nouv;
        }
        dist = fmax(dist,dist_nouv);
    }
    // On stock la distance avec la solution exacte au bout du vecteur u
    u[2*nbr_S+1]=dist;
    return u;
};

std::vector<double> BlackScholes::next_euler_expl(std::vector<double> u, float a,float b,float c, int taille, float pasT){
    std::vector<double> tab(taille,0.0);
    tab[0]=0;
    tab[taille-1]=u[taille-1];
    for (int j=1;j<(taille-1);++j){
        tab[j] = (b*pasT+1)*u[j]+c*pasT*u[j+1]+a*pasT*u[j-1];
    }
    return tab;
}

std::vector<double> BlackScholes::next_euler_impl(std::vector<double> u, float a,float b,float c, int taille, float pasT){
    float val_a = -a*pasT; // gauche
    float val_c = -c*pasT; // milieu
    float val_b = 1 -b*pasT; // droite

    // calcul de la solution sans avoir besoin d'inverser A
    std::vector<double> c_star(taille, 0.0);
    std::vector<double> d_star(taille, 0.0);
    
    c_star[0] = val_c/val_b;
    d_star[0] = u[0]/val_b;
    
    for (int i=1; i<taille; i++) {
        double m = 1.0 / (val_b - val_a * c_star[i-1]);
        c_star[i] = val_c * m;
        d_star[i] = (u[i] - val_a * d_star[i-1]) * m;
    }
    //std::vector<double> f(taille,0.0);
    for (int i=taille-2; i-- > 0; ) {
        u[i] = d_star[i] - c_star[i] * u[i+1];
    }
    return u;
}

std::vector<double> BlackScholes::next_Crank_Nicholson(std::vector<double> u, float a,float b,float c, int taille, float pasT, float theta){
    
    float val_a = -a*pasT*theta;// gauche
    float val_b = 1 -b*pasT*theta; // milieu
    float val_c = -c*pasT*theta; // droite
    
    std::vector<double> u_bis(taille, 0.0);
    // modification de u pour correspondre au membre de droite
    for (int j=1;j<(taille-1);++j){
        u_bis[j] = (b*(1-theta)*pasT+1)*u[j]+c*(1-theta)*pasT*u[j+1]+a*(1-theta)*pasT*u[j-1];
    }
    u_bis[0] = (b*(1-theta)*pasT+1)*u[0] + c*(1-theta)*pasT*u[1];
    u_bis[taille-1] = (b*(1-theta)*pasT+1)*u[taille-1] + b*(1-theta)*pasT*u[taille-2];
    
    // calcul de la solution sans avoir besoin d'inverser A
    std::vector<double> c_star(taille, 0.0);
    std::vector<double> d_star(taille, 0.0);
    
    c_star[0] = val_c/val_b;
    d_star[0] = u_bis[0]/val_b;
    
    for (int i=1; i<taille; i++) {
        double m = 1.0 / (val_b - val_a * c_star[i-1]);
        c_star[i] = val_c * m;
        d_star[i] = (u_bis[i] - val_a * d_star[i-1]) * m;
    }
    for (int i=taille-2; i-- > 0; ) {
        u_bis[i] = d_star[i] - c_star[i] * u_bis[i+1];
    }
    return u_bis;
}

void BlackScholes::write_elem_finie(float maxTemps,float maxS, int nbr_T, int nbr_S,int val){
    std::ofstream outputFile;
    outputFile.open("elem_finie_" + std::to_string(val)+"_"+ std::to_string(nbr_T)+".csv");
    std::vector<double> u = element_finie(maxTemps,maxS, nbr_T,nbr_S,val);
    for(int a=0;a<(2*nbr_S+1);++a){
        outputFile << u[a] << '\n' ;
    }
    outputFile.close();
};

void BlackScholes::write_exact_BS(float nbr_S, float max_S){
    std::ofstream outputFile;
    outputFile.open("exact_BS.csv");
    outputFile << "S"<< '\n' ;
    for(int j=0;j<(2*nbr_S+1);++j){
        float new_S0 = exp((j-nbr_S)*max_S/(nbr_S));
        change_S0(new_S0);
        float C_S = get_C();
        outputFile << C_S << '\n' ;
    }
    outputFile.close();
};

std::vector<double> BlackScholes::val_exact_BS(float nbr_S, float max_S){
    std::vector<double> u(2*nbr_S+1,0);
    float C_S =0;
    for(int j=0;j<(2*nbr_S+1);++j){
        float new_S0 = exp((j-nbr_S)*max_S/(nbr_S));
        change_S0(new_S0);
        C_S = get_C();
        u[j]=C_S;
    }
    return u;
};

float BlackScholes::get_d1_S_T(float S_0, float time){
    return (log(S_0/K)+(r+pow(sigma,2)/2)*time)/(sigma*sqrt(time));
};
float BlackScholes::get_d2_S_T(float S_0, float time){
    return get_d1_S_T(S_0, time)-sigma*sqrt(time);
};

std::vector<double> BlackScholes::val_exact_BS_T(float nbr_S, float max_S, float time_res){
    std::vector<double> u(2*nbr_S+1,0);
    for(int j=0;j<(2*nbr_S+1);++j){
        float new_S0 = exp((j-nbr_S)*max_S/(nbr_S));
        float x = new_S0*phi(get_d1_S_T(new_S0,time_res));
        float y = K*exp(-1*r*time_res)*phi(get_d2_S_T(new_S0,time_res));
        u[j]=x-y;
    }
    return u;
};

void BlackScholes::write_EDP_M(float maxTemps,float maxS, int nbr_T_0, int pas_T,int nbr_S,int val, int nbr_val){
    std::ofstream outputFile;
    outputFile.open("EDP_M_" + std::to_string(val)+".csv");
    float dist=0;
    for(int a = 0 ; a < nbr_val ; ++a){
        std::vector<double>u = element_finie(maxTemps,maxS, nbr_T_0+pas_T*a ,nbr_S,val);
        dist = u[2*nbr_S+1];
        outputFile << sqrt(dist) << '\n' ;
    }
    outputFile.close();
    std::cout << "fin EDP_M" << '\n';
};


void BlackScholes::write_EDP_N(float maxTemps,float maxS, int nbr_T, int pas_S, int nbr_S_0,int val, int nbr_val){
    std::ofstream outputFile;
    outputFile.open("EDP_N_" + std::to_string(val)+".csv");
    float dist =0;
    for(int a = 0 ; a < nbr_val ; ++a){
        std::vector<double> u = element_finie(maxTemps,maxS, nbr_T,nbr_S_0+pas_S*a ,val);
        dist = u[2*(nbr_S_0+pas_S*a)+1];
        outputFile << sqrt(dist) << '\n' ;
    }
    outputFile.close();
    std::cout << "fin EDP_N" << '\n';
};

std::vector<double> BlackScholes::next_American_impl(std::vector<double> u, float a,float b,float c, int taille, float pasT,float h, float L){
    float val_a = -a*pasT; // gauche
    float val_c = -c*pasT; // milieu
    float val_b = 1 -b*pasT; // droite
    
    
    // calcul de la solution sans avoir besoin d'inverser A
    std::vector<double> c_star(taille, 0.0);
    std::vector<double> d_star(taille, 0.0);
    
    c_star[0] = val_c/val_b;
    d_star[0] = u[0]/val_b;
    
    for (int i=1; i<taille; i++) {
        double m = 1.0 / (val_b - val_a * c_star[i-1]);
        c_star[i] = val_c * m;
        d_star[i] = (u[i] - val_a * d_star[i-1]) * m;
    }
    for (int i=taille-2; i-- > 0; ) {
        u[i] = d_star[i] - c_star[i] * u[i+1];
        double S = exp(-L+i*h);
        if(u[i]<K-S){
            u[i] = K-S ;
        }
    }
    return u;
}


std::vector<double> BlackScholes::next_Bermudan_impl(std::vector<double> u, float a,float b,float c, int taille, float pasT, float h, float L, int nStepsPerDate){
    float val_a = -a*pasT; // gauche
    float val_c = -c*pasT; // milieu
    float val_b = 1 -b*pasT; // droite
    
    
    // calcul de la solution sans avoir besoin d'inverser A
    std::vector<double> c_star(taille, 0.0);
    std::vector<double> d_star(taille, 0.0);
    
    c_star[0] = val_c/val_b;
    d_star[0] = u[0]/val_b;
    
    for (int i=1; i<taille; i++) {
        double m = 1.0 / (val_b - val_a * c_star[i-1]);
        c_star[i] = val_c * m;
        d_star[i] = (u[i] - val_a * d_star[i-1]) * m;
    }

    for (int i=taille-2; i-- > 0; ) {
        u[i] = d_star[i] - c_star[i] * u[i+1];
        // Compare to immediate exercise when allowed
        if( i % nStepsPerDate == 0){
            double S = exp(-L+i*h) ;
            if(u[i]<K-S){
                u[i] = K-S;
            }
        }
    }
    return u;
}
