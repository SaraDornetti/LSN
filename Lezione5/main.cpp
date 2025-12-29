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


double r (vector<double> x) {
   return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}


double p_uno_s (vector<double> x) { //calcolo direttamente il modulo quadro di psi
   return exp(-2*r(x))/M_PI;
}


double p_due_p (vector<double> x) {
   double rr = r(x);
   if (rr == 0.0) return 0.; // protezione contro divisione per zero
   return pow(rr*(x[2]/rr),2)*exp(-rr)/(32*M_PI);
}


void MCstep (vector<double>& x, double passo, Random& rnd, int& contatore, int config, int distr_type) {

   vector <double> dx = {0.,0.,0.};
   double alpha;

   if (distr_type==0){ //distr_type == 0 -> uniforme
      for (int i=0; i<3; i++) dx[i] = rnd.Rannyu(-passo, passo);
   }
   if (distr_type==1){ //distr_type == 0 -> Gaussiana
      for (int i=0; i<3; i++) dx[i] = rnd.Gauss(0.,passo);
   }

   // nuova configurazione candidta
   vector<double> newx(3);
   for (int i=0;i<3;i++) newx[i] = x[i] + dx[i];

   if (config==0) alpha = min(1., p_uno_s(newx)/p_uno_s(x)); //config == 0 -> 1s
   if (config==1) alpha = min(1., p_due_p(newx)/p_due_p(x)); //config == 1 -> 2p

   double acc = rnd.Rannyu();     
   if (acc<=alpha) {
      x = newx;
      contatore++;
   } 

   return;
}


void Equilibrazione (int n_equil, vector<double>& x, double passo, Random& rnd, int& contatore, int config, int distr_type, const string &filename, bool salva=false) {

   ofstream out;
   if (salva) {
       out.open(filename);
       if (!out) { cerr << "Errore apertura file " << filename << endl; exit(1); }
   }

   for (int i=0; i<n_equil; i++) {
      MCstep (x, passo, rnd, contatore, config, distr_type);
      if (salva) {
            out << i+1 << "\t" << r(x) << "\t" << x[0] << "\t" << x[1] << "\t" << x[2] << "\n";
        }
      }
   if (salva) out.close();
   return;

}



void BlockAverage (int M, int N, vector<double>& x, double passo, Random& rnd, int& contatore, int config, int distr_type, const string &filename) {
   int L = M/N;  // numero di lanci per blocco
   vector<double> ave(N,0.), ave2(N,0.), sum_prog(N,0.), sum2_prog(N,0.), err_prog(N,0.); // per media uniforme
  
   for (int i=0; i<N; i++) {
      double sum1=0;
      for (int j=0; j<L; j++) {
         MCstep (x, passo, rnd, contatore, config, distr_type);
         sum1 += r(x); 
      }
      ave[i] = sum1/L;           // media
      ave2[i] = ave[i]*ave[i];   // media al quadrato
  }   

   // calcolo medie progressive + errori
   for (int i=0; i<N; i++) {
      for (int j=0; j<=i; j++) {
          sum_prog[i]+=ave[j];         // somme medie progressive
          sum2_prog[i]+=ave2[j];       // somme medie al quadrato progressive
      }
    sum_prog[i]/=(i+1);
    sum2_prog[i]/=(i+1); 
    err_prog[i]=error(sum_prog, sum2_prog, i);
  }

  // scrivo su file blocco || media || errore
  ofstream out(filename);
  if (!out) { cerr << "Errore apertura file " << filename << endl; exit(1); }
  for (int i=0; i<N; i++) out << (i+1) << "\t" << sum_prog[i] << "\t" << err_prog[i] << "\n";
  out.close(); 
  return;
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


   int M=1000000;
   int N=100;
   

   // ===========================================
   // Configurazione 1s, distribuzione uniforme
   // ===========================================


   vector<double> x = {10., -8., 4.};
   double passo = 1.22;
   int contatore = 0;
   int n_equil = 10000;
   
   // termalizzazione
   Equilibrazione (n_equil, x, passo, rnd, contatore, 0, 0, "terma_U1.txt", true );
   cout << "Acceptance fraction during equil (U1): " << contatore/double(n_equil) << endl;
   // misura con blocking
   contatore = 0;
   BlockAverage(M, N, x, passo, rnd, contatore, 0, 0, "output_U1.txt");
   cout << "Acceptance fraction total (U1): " << contatore/double(M) << endl;




   // ===========================================
   // Configurazione 2p, distribuzione uniforme
   // ===========================================

   x = {-12., 23., -5.};
   passo = 2.8;
   contatore = 0;
   n_equil = 10000;

   // termalizzazione
   Equilibrazione (n_equil, x, passo, rnd, contatore, 1, 0, "terma_U2.txt", true );
   cout << "Acceptance fraction during equil (U2): " << contatore/double(n_equil) << endl;
   // misura con blocking
   contatore = 0;
   BlockAverage(M, N, x, passo, rnd, contatore, 1, 0, "output_U2.txt");
   cout << "Acceptance fraction total (U2): " << contatore/double(M) << endl;




   // ===========================================
   // Configurazione 1s, distribuzione Gaussiana
   // ===========================================

   x = {10., -8., 4.};;
   passo = 0.69;
   contatore = 0;
   n_equil = 10000;

   // termalizzazione
   Equilibrazione (n_equil, x, passo, rnd, contatore, 0, 1, "terma_G1.txt", true );
   cout << "Acceptance fraction during equil (G1): " << contatore/double(n_equil) << endl;
   // misura con blocking
   contatore = 0;
   BlockAverage(M, N, x, passo, rnd, contatore, 0, 1, "output_G1.txt");
   cout << "Acceptance fraction total (G1): " << contatore/double(M) << endl;


   // ===========================================
   // Configurazione 2p, distribuzione Gaussiana
   // ===========================================

   x = {-12., 23., -5.};
   passo = 1.75;
   contatore = 0;
   n_equil = 10000;

   // termalizzazione
   Equilibrazione (n_equil, x, passo, rnd, contatore, 1, 1, "terma_G2.txt", true );
   cout << "Acceptance fraction during equil (G2): " << contatore/double(n_equil) << endl;
   // misura con blocking
   contatore = 0;
   BlockAverage(M, N, x, passo, rnd, contatore, 1, 1, "output_G2.txt");
   cout << "Acceptance fraction total (G2): " << contatore/double(M) << endl;


rnd.SaveSeed();

return 0;

} 





