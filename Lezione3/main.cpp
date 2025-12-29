#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "random.h"

using namespace std;

double error(const vector<double>& AV, const vector<double>& AV2, int n) {
    if (n == 0) {
        return 0.0;
    } else {
        return sqrt((AV2[n] - AV[n] * AV[n]) / n);
    }
}

double N(double x) {
    return 0.5 * (1. + erf(x / sqrt(2.)));
}

double call_black_scholes(double S0, double K, double T, double r, double sigma) {
    double d1 = 1./(sigma * sqrt(T)) * (log(S0 / K) + (r + (pow(sigma,2)) / 2.) * T);
    double d2 = d1 - sigma * sqrt(T);
    double C = S0 * N(d1) - K * exp(-r * T) * N(d2);
   return C;
}

double put_black_scholes(double S0, double K, double T, double r, double sigma) {
   double d1 = 1./(sigma * sqrt(T)) * (log(S0 / K) + (r + (pow(sigma,2)) / 2.) * T);
   double d2 = d1 - sigma * sqrt(T);
   double P = S0 *(N(d1) - 1.) - K * exp(-r * T) * (N(d2)-1.);
  return P;
}

double S (double S0, double mu, double sigma, double t, double Z ) {
   return S0 * exp((mu-(0.5 * sigma * sigma)) * t + sigma * Z * sqrt(t));
}


//io in questo caso uso intervalli della stessa lunghezza, in realtà basta che i tempi siano in ordine
double S_recursive (double S0, double mu, double sigma, double t, Random &rnd) {

   int n_intervals = 100; //numero di intervalli
   double t_interval = t/n_intervals;
   double S = S0;

   for (int i=0; i<n_intervals; i++) {
      double z = rnd.Gauss(0,1);
      S = S * exp((mu-(0.5 * sigma * sigma)) * (t_interval*(i+1)-t_interval*i) + sigma * z * sqrt((t_interval*(i+1)-t_interval*i)));
   }

   return S;

}

double C (double r, double T, double S, double K) {
  return exp(-r*T)*max(0., S - K);
}


double P (double r, double T, double S, double K) {
   return exp(-r*T)*max(0., K - S);
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

   
// data
double S0 = 100.; //costo al tempo zero
double K = 100.; //a quanto lo vendo: al prezzo al tempo zero
double T = 1.; //tempo
double r = 0.1; //di quanto si alza il prezzo
double sigma = 0.25; //volatilita'
// double t=0;
   


// data blocking
int M=10000;
int N=100;
int L=M/N;


/******************************************************************************
Sampling direct GBM
*******************************************************************************/

vector<double> ave_P, ave_C, ave2_P, ave2_C, sum_prog_P, sum2_prog_P, err_prog_P, sum_prog_C, sum2_prog_C, err_prog_C ;


for (int i=0; i<N; i++) {

   double sumP=0;
   double sumC=0;

   for (int j=0; j<L; j++) {
   double Z = rnd.Gauss(0,1);
   double St = S (S0, r, sigma, T, Z ); //St e' il prezzo al tempo t
   sumC += C(r, T, St, K); 
   sumP += P(r, T, St, K);
   }

   ave_C.push_back(sumC/L); 
   ave2_C.push_back(ave_C[i]*ave_C[i]);

   ave_P.push_back(sumP/L); 
   ave2_P.push_back(ave_P[i]*ave_P[i]); 

}



for (int i=0; i<N; i++) {
   double sum2_C=0;
   double sum3_C=0;

   double sum2_P=0;
   double sum3_P=0;

   for (int j=0; j<i+1; j++) {
      sum2_C += ave_C[j];
      sum3_C += ave2_C[j];

      sum2_P += ave_P[j];
      sum3_P += ave2_P[j];
   }

   sum_prog_C.push_back(sum2_C/(i+1));
   sum2_prog_C.push_back(sum3_C/(i+1));
   err_prog_C.push_back(error(sum_prog_C, sum2_prog_C, i));

   sum_prog_P.push_back(sum2_P/(i+1));
   sum2_prog_P.push_back(sum3_P/(i+1));
   err_prog_P.push_back(error(sum_prog_P, sum2_prog_P, i));
   
}



ofstream outputFile("outputC.txt");  // Crea o sovrascrive il file
ofstream outputFile2("outputP.txt");  // Crea o sovrascrive il fi

if (!outputFile) {
   cerr << "Errore nell'apertura del file!" << endl;
   return 1;  // Esce con errore
}
if (!outputFile2) {
   cerr << "Errore nell'apertura del file!" << endl;
   return 1;  // Esce con errore
}


for (int i=0; i<N; i++){
   outputFile << (i+1) << "\t" << sum_prog_C[i] << "\t" << err_prog_C[i] << endl;
   outputFile2 << (i+1) << "\t" << sum_prog_P[i] << "\t" << err_prog_P[i] << endl;
}

outputFile.close();  // Chiudi il file
outputFile2.close();  // Chiudi il file
    





/******************************************************************************
Sampling discretized GBM
*******************************************************************************/

vector<double> ave_P2, ave_C2, ave2_P2, ave2_C2, sum_prog_P2, sum2_prog_P2, err_prog_P2, sum_prog_C2, sum2_prog_C2, err_prog_C2 ;


for (int i=0; i<N; i++) {

   double sumP=0;
   double sumC=0;

   for (int j=0; j<L; j++) {
   double St = S_recursive (S0, r, sigma, T, rnd); //St e' il prezzo al tempo t
   sumC += C(r, T, St, K); 
   sumP += P(r, T, St, K);
   }

   ave_C2.push_back(sumC/L); 
   ave2_C2.push_back(ave_C2[i]*ave_C2[i]);

   ave_P2.push_back(sumP/L); 
   ave2_P2.push_back(ave_P2[i]*ave_P2[i]); 

}



for (int i=0; i<N; i++) {
   double sum2_C=0;
   double sum3_C=0;

   double sum2_P=0;
   double sum3_P=0;

   for (int j=0; j<i+1; j++) {
      sum2_C += ave_C2[j];
      sum3_C += ave2_C2[j];

      sum2_P += ave_P2[j];
      sum3_P += ave2_P2[j];
   }

   sum_prog_C2.push_back(sum2_C/(i+1));
   sum2_prog_C2.push_back(sum3_C/(i+1));
   err_prog_C2.push_back(error(sum_prog_C2, sum2_prog_C2, i));

   sum_prog_P2.push_back(sum2_P/(i+1));
   sum2_prog_P2.push_back(sum3_P/(i+1));
   err_prog_P2.push_back(error(sum_prog_P2, sum2_prog_P2, i));
   
}



ofstream outputFile3("outputC2.txt");  // Crea o sovrascrive il file
ofstream outputFile4("outputP2.txt");  // Crea o sovrascrive il fi

if (!outputFile3) {
   cerr << "Errore nell'apertura del file!" << endl;
   return 1;  // Esce con errore
}
if (!outputFile4) {
   cerr << "Errore nell'apertura del file!" << endl;
   return 1;  // Esce con errore
}


for (int i=0; i<N; i++){
   outputFile3 << (i+1) << "\t" << sum_prog_C2[i] << "\t" << err_prog_C2[i] << endl;
   outputFile4 << (i+1) << "\t" << sum_prog_P2[i] << "\t" << err_prog_P2[i] << endl;
}

outputFile3.close();  // Chiudi il file
outputFile4.close();  // Chiudi il file


rnd.SaveSeed();

return 0;

} 




