#ifndef BlackScholes_hpp
#define BlackScholes_hpp

#include <stdio.h>
#include <fstream>

class BlackScholes
{
    
// Méthodes
public:
    
    BlackScholes(float r, float s, float T, float S0, float K );
    
    
    // Pricing d'un Call Européen Black Scholes par la méthode de Monte Carlo
        //Calcul de l'espérance du payoff
    float get_d1();
    float get_d2();
    float get_C();
    float phi(float d);
    float get_Z();
    float get_X();
    void change_S0(float s);
    
        // Evaluation des performances de la méthode de Monte Carlo
    float* MC(int n);
    float eval_time_MC(int n, int m);
    float* eval_val_MC(int n, int m);
    void write_MC(int m, int N, int k);
    
    
    // EDP - Éléments finies
    
    // algorithme principale
    std::vector<double> element_finie(float maxTemps,float maxS, int nbr_T, int nbr_S,int val);
    // differentes méthodes étudiées
    std::vector<double> next_euler_impl(std::vector<double> u, float a,float b,float c, int taille, float pasT);
    std::vector<double> next_euler_expl(std::vector<double> u, float a,float b,float c, int taille, float pasT);
    std::vector<double> next_Crank_Nicholson(std::vector<double> u, float a,float b,float c, int taille, float pasT, float theta);
    // ecriture du fichier, la dernière valeur est la norme
    void write_elem_finie(float maxTemps,float maxS, int nbr_T, int nbr_S, int val);
    
    
    void write_exact_BS(float nbr_S, float max_S);
    std::vector<double> val_exact_BS(float nbr_S, float max_S);
    
    void write_EDP_M(float maxTemps,float maxS, int nbr_T_0, int pas_T,int nbr_S,int val, int nbr_val);
    void write_EDP_N(float maxTemps,float maxS, int nbr_T, int pas_S, int nbr_S_0,int val, int nbr_val);
    std::vector<double> val_exact_BS_T(float nbr_S, float max_S, float time_res);
    
    
    // American Option et Bermudan
    std::vector<double> next_American_impl(std::vector<double> u, float a,float b,float c, int taille, float pasT,float h, float L);
    std::vector<double> next_Bermudan_impl(std::vector<double> u, float a,float b,float c, int taille, float pasT,float h, float L, int nStepsPerDate);
    void write_elem_finie_TP6(float maxTemps,float maxS, int nbr_T, int nbr_S,int val);
    float get_d1_S_T(float S_0, float time);
    float get_d2_S_T(float S_0, float time);
    
// Attributs du Modèle
private:
    float r;
    float sigma;
    float T;
    float S0;
    float K;
    
};


#endif /* BlackScholes_hpp */
