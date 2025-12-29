/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include "system.h"

using namespace std;

int main (int argc, char *argv[]){

  // CREO E INIZIALIZZO IL SISTEMA

  System SYS;
  SYS.initialize();
  SYS.initialize_properties();
  SYS.block_reset(0);


  // EQUILIBRAZIONE

  
  /*for(int i=0; i < SYS.get_nbl(); i++){ //loop over blocks
    for(int j=0; j < SYS.get_nsteps(); j++){ //loop over steps in a block
      SYS.step();
      SYS.measure();
    }
    SYS.averages(i+1);
    SYS.block_reset(i+1);
  }
  SYS.finalize();*/
  

  // SIMULAZIONE

  int N_eq = 6000; 

  int n_stepT = 20;
  double T = 2.0;
  double T_min = 0.5;
  double dT = (T-T_min)/n_stepT;

  for (int n=0; n<=n_stepT; n++){

    SYS.set_temp(T);
    SYS.reset_averages();

    //equilibrazione
    for (int m = 0; m < N_eq; m++) {
      SYS.step();
    }

    for(int i=0; i < SYS.get_nbl(); i++){ //loop over blocks
      SYS.block_reset(i+1); 
      for(int j=0; j < SYS.get_nsteps(); j++){ //loop over steps in a block
        SYS.step();
        SYS.measure();
      }
      SYS.averages(i+1);
      SYS.block_reset(i+1);
    }
    SYS.finalize();

    T -= dT;

  }

  return 0;

}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
