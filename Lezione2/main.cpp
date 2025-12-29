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

int M=10000;
int N=100;
int L=M/N;

vector<double> ave, ave2, sum_prog, sum2_prog, err_prog;

//calcolo integrale con distribuzione uniforme

for (int i=0; i<N; i++) {
    double sum1=0;
    for (int j=0; j<L; j++) {
        sum1 += (M_PI/2)*cos(M_PI*rnd.Rannyu()/2); //controlla che venga double
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

ofstream outputFile("output.txt");  // Crea o sovrascrive il file

   if (!outputFile) {
      cerr << "Errore nell'apertura del file!" << endl;
      return 1;  // Esce con errore
   }

   for (int i=0; i<N; i++){
      outputFile << i+1 << "\t" << sum_prog[i] << "\t" << err_prog[i] << endl;
}

outputFile.close();  // Chiudi il file
    


// Calcolo integrale con distribuzione diversa --> importance sampling

ofstream outputFile_test("output_test.txt");
for (int i=0; i<100000; i++) {
   outputFile_test << i+1 << "\t" << rnd.Line() << endl;
}
outputFile_test.close();


vector<double> ave_is, ave2_is, sum_prog_is, sum2_prog_is, err_prog_is;

for (int i=0; i<N; i++) {
   double sum1=0;
   for (int j=0; j<L; j++) {
      double x = rnd.Line();
       sum1 += (M_PI/2)*cos(M_PI*x/2)/(-2*x+2); //controlla che venga double
   }
   ave_is.push_back(sum1/L); 
   ave2_is.push_back(ave_is[i]*ave_is[i]);
}

for (int i=0; i<N; i++) {
   double sum2=0;
   double sum3=0;
   for (int j=0; j<i+1; j++) {
      sum2 += ave_is[j];
      sum3 += ave2_is[j];
   }

   sum_prog_is.push_back(sum2/(i+1));
   sum2_prog_is.push_back(sum3/(i+1));
   err_prog_is.push_back(error(sum_prog_is, sum2_prog_is, i));
   
}

ofstream outputFile_is("output_is.txt");  // Crea o sovrascrive il file

   if (!outputFile_is) {
      cerr << "Errore nell'apertura del file!" << endl;
      return 1;  // Esce con errore
   }

   for (int i=0; i<N; i++){
      outputFile_is << i+1 << "\t" << sum_prog_is[i] << "\t" << err_prog_is[i] << endl;
}

outputFile_is.close();  // Chiudi il file

rnd.SaveSeed();

return 0;

}

