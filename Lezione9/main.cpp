#include <iostream>
#include "system.h"

using namespace std;

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


  System SYS;
  SYS.initialize(rnd);

  for(int i=0; i < SYS.get_ngenerations(); i++){ 
    SYS.new_generation(rnd);
  }

  SYS.finalize();

  return 0;
}

