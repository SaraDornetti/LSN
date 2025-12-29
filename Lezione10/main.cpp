#include <iostream>
#include <mpi.h>
#include "system.h"

using namespace std;

int main (int argc, char *argv[]){

   // MPI initialization
   int size, rank; 
   MPI_Init(&argc, &argv); 
   MPI_Comm_size(MPI_COMM_WORLD, &size); 
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   
   bool enable_migration = false; 
   int N_migr = 50;  //ogni N_migr generazioni si effettuano le migrazioni 

  // Random number generator initialization
  Random rnd;
  string seed_file = "seed.in";
  string primes_file = "Primes";
  rnd.Initialize(seed_file, primes_file, rank); 

  System SYS;
  SYS.initialize(rnd);
  
  int n_gen = SYS.get_ngenerations();

  for(int i_gen=0; i_gen < n_gen; i_gen++){ 
    SYS.new_generation(rnd, rank);
    if (enable_migration && i_gen % N_migr == 0 && i_gen > 0) {
      SYS.migration(rank, size);
      //SYS.selection(); // l'ho chiamata selection anche se in realtà riordina e basta
  }
}

Individuo best_individual = SYS.get_individuo(0);
best_individual.compute_fitness();
double local_best_distance = best_individual.get_fitness();

// This vector will contain the best distances for each process
vector<double> all_best_distances(size);
// All processes send their local_best_distance to rank 0, which collects them in all_best_distances
MPI_Gather(&local_best_distance, 1, MPI_DOUBLE,
           all_best_distances.data(), 1, MPI_DOUBLE,
           0, MPI_COMM_WORLD);

// Search the rank corresponding to the best path
int best_rank = 0;
if (rank == 0) {
    double min_dist = all_best_distances[0];
    for (int i = 1; i < size; i++) {
        if (all_best_distances[i] < min_dist) {
            min_dist = all_best_distances[i];
            best_rank = i;
        }
    }
}

// All processes must know which is the rank with the best path
MPI_Bcast(&best_rank, 1, MPI_INT, 0, MPI_COMM_WORLD);

// The best rank saves his results
if (rank == best_rank) {
   SYS.write_configuration();
}

MPI_Finalize();
  
return 0;

}

