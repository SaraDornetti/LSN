#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include "random.h"

using namespace std;

double error(const vector<double>& AV, const vector<double>& AV2, int n) {
    if (n == 0) {
        return 0.0;
    } else {
        return sqrt((AV2[n] - AV[n] * AV[n]) / n);
    }
}

double psi_quadro (double x, double sigma, double mu) { //calcolo direttamente il modulo quadro di psi
   return pow(exp(-pow(x - mu, 2) / (2*pow(sigma, 2))) + exp(-pow(x + mu, 2) / (2*pow(sigma, 2))),2);
}

double wf(double x, double sigma, double mu){
   return exp(-pow(x - mu, 2) / (2*pow(sigma, 2))) + exp(-pow(x + mu, 2) / (2*pow(sigma, 2)));
}

double laplace(double x, double sigma, double mu){
   double factor = -1./(sigma*sigma);
   double frac_minus = pow((x + mu)/sigma,2);
   double frac_plus = pow((x - mu)/sigma,2);
   double exp_minus= exp(-1. * frac_minus/2);
   double exp_plus = exp(-1. * frac_plus/2);
   return factor * ((1. - frac_minus) * exp_minus + (1. - frac_plus) * exp_plus);
}

double kinetic_energy(double x, double sigma, double mu){
   return - 1./2. * laplace(x, sigma, mu)/wf(x, sigma, mu); // Natural units
}

double potential_energy(double x){
   return pow(x,4)-2.5*pow(x,2); // Natural units
}

double Metropolis_aveEnergy(double mu, double sigma, Random &rnd) {

   int M=100000; // lanci del Metropolis_aveEnergy

   double ave_energy=0.;
   double x=0.;
   double passo = 2.;

   for (int i=0; i<M; i++) {

      double x1=rnd.Rannyu(-passo, passo);
      double alpha = min(1., psi_quadro(x+x1, sigma, mu)/psi_quadro(x, sigma, mu));
      double acc = rnd.Rannyu();
      
      if (acc<=alpha) x=x+x1;
        
      ave_energy += kinetic_energy(x, sigma, mu)+potential_energy(x); // ho già normalizzato nella funzione E_kin
      
   }   
   
   return ave_energy/M;

}


double Acceptance_Ratio_SA (double energy_old, double energy_new, double T) {
   return min(1., exp(-1./T*(energy_new-energy_old)));
}

int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

   
   
   double mu_start=1.;
   double sigma_start=1.;


   double energy_old = Metropolis_aveEnergy(mu_start, sigma_start, rnd);
   double energy_new = 0.;

   int n_T = 200;
   int n_cicli = 20;
   double step = 0.5;

   ofstream outputFile_SA("output_SA.txt");  // Crea o sovrascrive il file
   
   if (!outputFile_SA) {
      cerr << "Errore nell'apertura del file!" << endl;
      return 1;  // Esce con errore
   } 


   for (int i=0; i<n_T; i++) {

      cout << "STEP " << i+1 << "\t" << "dim_step = " << step << endl;

      int conteggio = 0;

      double T = 3./(1+4*i);

      for (int k=0; k<n_cicli; k++) { 

         double sigma = sigma_start + rnd.Rannyu(-step, step);
         double mu = mu_start + rnd.Rannyu(-step, step);

         energy_new = Metropolis_aveEnergy(mu, sigma, rnd);

         if (rnd.Rannyu()<Acceptance_Ratio_SA(energy_old, energy_new, T)) {
            sigma_start=sigma;
            mu_start=mu;
            energy_old=energy_new;
            conteggio ++;
            cout << "mu = " << mu << "\t" << "sigma = " << sigma << "\t" << "energia = " << energy_old << endl;
            outputFile_SA << T << "\t" << mu_start << "\t" << sigma_start << "\t" << energy_old << endl;
         }

         
      
      }

      cout << "Per la temperatura T = " << T << "\t" << "accettanza = " << double(conteggio)/n_cicli << "%" << endl << endl;
      
      if (double(conteggio)/n_cicli<0.40 && step > 0.01) {
         step = step*0.9;
      }
      if (double(conteggio)/n_cicli>0.60 && step < 2.) {
         step = step*1.1;
      }


   }

   
rnd.SaveSeed();

return 0;

} 






