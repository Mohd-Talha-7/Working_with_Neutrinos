#include "Math/Integrator.h"
#define c 10

double func(double E_nu) 
{
    
   /* product of the fission fractions & neutrino flux per fission */
                                                                     
    double Fi_Si1 = 0.564 * exp(0.870 - 0.160*E_nu - 0.0910*pow(E_nu,2));  //  U235
    double Fi_Si2 = 0.304 * exp(0.896 - 0.239*E_nu - 0.0981*pow(E_nu,2));  //  Pu239
    double Fi_Si3 = 0.076 * exp(0.976 - 0.162*E_nu - 0.0790*pow(E_nu,2));  //  U238
    double Fi_Si4 = 0.056 * exp(0.793 - 0.080*E_nu - 0.1085*pow(E_nu,2));  //  Pu241
    
    /*Oscillation parameters Normal hierarchy Nu-fit 2021 */
    
    double dm_12 = 7.42e-5; /*eV2*/ 
    double dm_13 = 2.515e-3; /*eV2*//*check Nufit*/
    double dm_23 = dm_13-dm_12;/*eV2*/
    double s_12 = 0.304;
    double s_23 = 0.573;
    double s_13 = 0.02220;
    double th12= 33.44 * 4.0*atan(1.0)/180;
    double th23= 49.20 * 4.0*atan(1.0)/180;
    double th13= 8.57  * 4.0*atan(1.0)/180;
    double D21 = 1.27 * dm_12 * 52476/E_nu;
    double D31 = 1.27 * dm_13 * 52476/E_nu;
    double D32 = 1.27 * dm_23 * 52476/E_nu;
    
    
    double P_21 = pow(cos(th13),4) * pow(sin(2.0*th12),2) * pow(sin(D21),2);
    double P_31 = pow(cos(th12),2) * pow(sin(2.0*th13),2) * pow(sin(D31),2);
    double P_32 = pow(sin(th12),2) * pow(sin(2.0*th13),2) * pow(sin(D32),2);
    
    double P_ee = 1.0 - P_21 - P_31 - P_32;     // expression for survival probability ref[65]
    
    double cross = 0.0952*(E_nu-1.29)*sqrt(pow((E_nu-1.29),2)-pow(0.511,2))*pow(10,-9);  // cross section [65], the order of the number of protons has been absorbed in the 
                                                                                         // cross section expression and therefore 10^{-9}                           
    
     double fi_ei = 0.564 * 202.36  + 0.076 * 205.99 + 0.304 * 211.12 +  0.056 * 214.26; // fission fraction & thermal energy/fission
     
     double Unit_conv = 6241506479963.2* 1e9; // 1 joule = 6.2415 x 10^12 MeV, 1 GW  --> 1e9
     
     double L_r = (1.0/(4.0* TMath::Pi() *pow(5247600,2))); // Mean of the baseline = 5247600 cm^2, 
   
     double W_r = (35.8/fi_ei) * Unit_conv;  // Total thermal power = 35.8 GW 
    
    
    double f= L_r * W_r * (Fi_Si1 + Fi_Si2 + Fi_Si3 + Fi_Si4)* cross *6.75 * P_ee;    // event rate expression taking probability into account [main article, pg. 191]
    
    return f;
}

int integ4() {
    
   // Set default tolerances for all integrators.
    ROOT::Math::IntegratorOneDimOptions::SetDefaultAbsTolerance(1.E-6);
    ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1.E-6); 
   
   ROOT::Math::Functor1D wf(&func);
   //ROOT::Math::Integrator ig(ROOT::Math::IntegrationOneDim::kLEGENDRE);
   ROOT::Math::Integrator ig(ROOT::Math::IntegrationOneDim::kGAUSS);
   //ROOT::Math::Integrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
   ig.SetFunction(wf);
    
   double Emin = 1.8, Emax = 10.0; 
      
      double val = ig.Integral(Emin, Emax);  
      
      val = val * 60.0 * 60.0 * 24.0; // changing the rate --> /sec to /day
    
    cout << "integral result is " << val <<endl; // it is giving the result 354.937. It seems we want a correction factor of about 0.23 to get 83 IBD events.
    
    return 0;
}

