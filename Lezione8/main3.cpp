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
   return pow(exp(-(x-mu)*(x-mu)/(2*sigma*sigma))+exp(-(x+mu)*(x+mu)/(2*sigma*sigma)),2);
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

   int M=10000; // lanci del Metropolis_aveEnergy

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

   
   ifstream file("output_SA.txt");
    if (!file) {
        cerr << "Errore nell'apertura del file.\n";
        return 1;
    }

    double a, b, c, d;
    vector<double> mu_vec, sigma_vec;

    while (file >> a >> b >> c >> d) {
        mu_vec.push_back(b); 
        sigma_vec.push_back(c);
    }

    file.close();


   int M=10000;
   int N=100;
   int L=M/N;

   
   double x = 0.;
   double mu, sigma=0;

   int contatore = 0;
   
   ofstream outputFile("output_energie.txt");  // Crea o sovrascrive il file
   if (!outputFile) {
      cerr << "Errore nell'apertura del file!" << endl;
      return 1;  // Esce con errore
   }

   for (int k=0; k<mu_vec.size(); k++) {

      mu = mu_vec[k];
      sigma = sigma_vec[k];

      vector<double> ave, ave2, sum_prog, sum2_prog, err_prog;

      for (int i=0; i<N; i++) {

         double sum1=0;
         double passo = 2.2;

         for (int j=0; j<L; j++) {

            double x1=rnd.Rannyu(-passo, passo);
        
            double alpha = min(1., psi_quadro(x+x1, sigma, mu)/psi_quadro(x, sigma, mu));
            double acc = rnd.Rannyu();
      

            if (acc<=alpha) {
               x=x+x1;
               contatore ++;
            }

            sum1 += kinetic_energy(x, sigma, mu)+potential_energy(x); // ho già normalizzato nella funzione E_kin
       
         }

         ave.push_back(sum1/L); 
         ave2.push_back(ave[i]*ave[i]);

      }



      for (int i=0; i<N; i++) {

         double sum2=0;
         double sum3=0;
   
         for (int j=0; j<i+1; j++) {
   
            sum2 += ave[j];
            sum3 += ave2[j];
   
         }
      
         sum_prog.push_back(sum2/(i+1));
         sum2_prog.push_back(sum3/(i+1));
         err_prog.push_back(error(sum_prog, sum2_prog, i));
         
      }
      
      outputFile << k << "\t" << sum_prog[N-1] << "\t" << err_prog[N-1] << endl;   

   }
   
   outputFile.close();  // Chiudi il file
       
   rnd.SaveSeed();
   return 0;

} 


