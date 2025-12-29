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


   int M=100000;
   int N=100;
   int L=M/N;
   
   vector<double> ave, ave2, sum_prog, sum2_prog, err_prog;
   
   double x = 0.;
   double mu = 0.77296;
   double sigma = 0.627109;
   double passo = 2.4;

   int contatore = 0;

   // TERMALIZZAZIONE
   for(int i = 0; i < 5000; i++){
      double x1 = rnd.Rannyu(-passo, passo);
      double alpha = min(1., psi_quadro(x+x1,sigma,mu)/psi_quadro(x,sigma,mu));
      if(rnd.Rannyu() <= alpha) x += x1;
   }

   ofstream outputFile_wf("output_wf.txt");  // Crea o sovrascrive il file
   
   if (!outputFile_wf) {
      cerr << "Errore nell'apertura del file!" << endl;
      return 1;  // Esce con errore
   } 

   for (int i=0; i<N; i++) {

       double sum1=0;

       for (int j=0; j<L; j++) {

         double x1=rnd.Rannyu(-passo, passo);
        
         double alpha = min(1., psi_quadro(x+x1, sigma, mu)/psi_quadro(x, sigma, mu));
         double acc = rnd.Rannyu();
      

         if (acc<=alpha) {
            x=x+x1;
            contatore ++;
         }

         outputFile_wf << i << "\t" << x << endl;
         sum1 += kinetic_energy(x, sigma, mu)+potential_energy(x); // ho già normalizzato nella funzione E_kin
       
      }

       ave.push_back(sum1/L); 
       ave2.push_back(ave[i]*ave[i]);
   }
   
   
   outputFile_wf.close();  // Chiudi il file

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
   

   ofstream outputFile("output_s.txt");  // Crea o sovrascrive il file
   
      if (!outputFile) {
         cerr << "Errore nell'apertura del file!" << endl;
         return 1;  // Esce con errore
      }
   
      for (int i=0; i<N; i++){
         outputFile << L*i << "\t" << sum_prog[i] << "\t" << err_prog[i] << endl;
   }
   
   outputFile.close();  // Chiudi il file
       
   
   cout << contatore/double(M) << endl;



rnd.SaveSeed();

return 0;

} 







