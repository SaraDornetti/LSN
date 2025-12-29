#ifndef __System__
#define __System__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <armadillo>
#include <stdlib.h> //exit
#include <mpi.h>
#include "individuo.h"
#include "random.h"

using namespace std;
using namespace arma;

class System {

private:
  int _city_shape;        // Type of city (e.g., circle, square,...)
  int _ncities;           // Number of cities
  int _nindividui;        // numero di individui in una generazione
  int _ngenerations;      // numero totale di generazioni
  double _crossover_rate; // rate con cui effettuare crossover (diminuisce aumentando il numero della generazione)
  double _mutation_rate;    // rate con cui effettuare mutazione (diminuisce aumentando il numero della generazione)
  double _p;                // esponente da mettere in selection  
  int _generation_index;    // indice della generazione
  field <Individuo> _individuo; // Field of particle objects representing the system
  field <Individuo> _individuo_figlio;
  mat _distance_matrix;

public: // Function declarations
  void initialize(Random& _rnd);          // Initialize system properties
  void initialize_circle(Random& _rnd);   // Inizializza un sistema con le città disposte in posizioni random su una circonferenza unitaria
  void initialize_square(Random& _rnd);   // Inizializza un sistema con le città disposte in posizioni random all'interno di un quadrato di lato uno
  void initialize_Italy();   // Inizializza un sistema con le città disposte in posizioni random all'interno di un quadrato di lato uno
  void new_generation(Random& _rnd, int rank);      // parte dalla generazione attuale e implementa crossover/mutazioni
  int get_ngenerations();
  void selection();           // Seleziona individui tra cui fare crossover
  void crossover(Random& _rnd, int ind1, int ind2); 
  void write_configuration(); // Write system configuration 
  void finalize();            // Finalize system and clean up
  
  void mean_fitness(int rank); //in cui calcolo la media dei fitness a ogni generazione e salvo in un file --> apri file in initialize()
  //usa generation index!!!!!!!!

  void migration(int rank, int size);
  Individuo get_individuo(int idx);


};

#endif // __System__
