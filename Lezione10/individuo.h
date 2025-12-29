#ifndef __Individuo__
#define __Individuo__

#include <armadillo>
#include "random.h"

using namespace std;
using namespace arma;

class Individuo {

private:
    int _ncities;           // numero di città
    mat _distance_matrix;
    double _fitness;
    //forse posso mettere uvec
    uvec _x;                 // Current position vector 
  

public: // Function declarations
    void initialize(Random& _rnd);                            // Initialize particle properties
    void check();                                 // controlla che la sequenza sia permessa
    void swap(int i, int j);                      // porta la città in i-esima posizione in j-esima 
    void setncities(int ncities);
    double get_fitness() const;                          // Get the fitness of the ind
    void compute_fitness();                       // compute the fitness of the ind
    uvec get_x() const;
    void set_x(const uvec &x);
    void set_distancematrix(const mat &distance_matrix);
    
    void pair_permutation(Random& _rnd);
    void shift(Random& _rnd);
    void contiguous_permutation(Random& _rnd);
    void inversion(Random& _rnd);
    
    
};


#endif // __Individuo__
 