#include <iostream>
#include <bits/stdc++.h>

using namespace std;

// Initial parameters definitions

// Diffusitivities
const double sig_ur = 2.8e-10;
const double sig_uz = 1.1e-9;
const double sig_vr = 2.32e-9;
const double sig_vz = 6.97e-9;

// Respiration kinetic parameters
const double T_ref = 293.15;
const double R_g = 8.314;
    // Oxygen consumption
double V_mu;
const double V_mu_ref = 2.39e-4;
const int Ea_vmu_ref = 80200;
    // Fermentative carbon dioxide production
double V_mfv;
const double V_mfv_ref = 1.61e-4;
const int Ea_vmfv_ref = 56700;

// Michaelis-Menten constants
const double K_mu = 0.4103;
const double K_mv = 27.2438;
const double K_mfu = 0.1149;

// Respiration constant
const double r_q = 0.97;

// Convective mass transfer coefficients
const double rho_u = 7e-7;
const double rho_v = 7.5e-7;

// Ambient conditions
const int p_atm = 101300;
double T;
double Cu_amb;
double Cv_amb;

// Parameters depend on storage conditions:
// storage = 1 if Orchard
// storage = 2 if Shelf life
// storage = 3 if Refrigerator
// storage = 4 if Precooling
// storage = 5 if Disorder inducing
// storage = 6 if Optimal CA

void setParameters(int storage){
    int T_cel;
    double eta_u;
    double eta_v;

    switch (storage) {
        case 1: T_cel = 25; eta_u = 20.8; eta_v = 0.04; break;
        case 2: T_cel = 20; eta_u = 20.8; eta_v = 0;    break;
        case 3: T_cel = 7;  eta_u = 20.8; eta_v = 0;    break;
        case 4: T_cel = -1; eta_u = 20.8; eta_v = 0;    break;
        case 5: T_cel = -1; eta_u = 2;    eta_v = 5;    break;
        case 6: T_cel = -1; eta_u = 2;    eta_v = 0.7;  break;
        default: cout << "Error! No valid storage conditions are specified"; break;
    }
    
    T = T_cel + 273.15;
    Cu_amb = (p_atm*eta_u/100)/(R_g*T);
    Cv_amb = (p_atm*eta_v/100)/(R_g*T);

    V_mu = V_mu_ref*exp(Ea_vmu_ref/R_g*(1/T_ref-1/T));
    V_mfv = V_mfv_ref*exp(Ea_vmfv_ref/R_g*(1/T_ref-1/T));

}