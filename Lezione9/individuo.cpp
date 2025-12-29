#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>
#include "individuo.h"

using namespace std;

void Individuo :: initialize(Random& _rnd){

   _x.set_size(_ncities); 
   for (int i=0; i<_ncities; i++) {
      _x(i)=i;
   }
   for (int i=0; i<_ncities; i++) {
      int j = static_cast<int>(_rnd.Rannyu(1,_ncities));
      int k = static_cast<int>(_rnd.Rannyu(1,_ncities));
      this->swap(j,k);
   }
   this->check();
   this->compute_fitness();
   return;
}


void Individuo::check() {
   // 1. Controlla che la prima città sia 0
   if (_x(0) != 0) {
       cout << "Errore: la prima città non è 0." << endl;
   }

   // 2. Controlla che non ci siano città ripetute
   vector<int> cities(_x.begin(), _x.end());
   sort(cities.begin(), cities.end());
   for (int i = 0; i < _ncities - 1; i++) {
       if (cities[i] == cities[i + 1]) {
           cout << "Errore: città ripetuta nel percorso." << endl;
           break;
       }
   }

   // 3. Controlla che siano presenti tutte le città da 0 a _ncities-1
   if (cities.front() != 0 || cities.back() != _ncities - 1) {
       cout << "Attenzione: il percorso non contiene tutte le città corrette." << endl;
   }
}



void Individuo :: swap(int i, int j) {
   int k = _x(i);
    _x(i) = _x(j);
    _x(j) = k;

    return;
}


void Individuo :: setncities(int ncities){
   _ncities = ncities;
   return;
}


double Individuo :: get_fitness() const{
   return _fitness;
}

void Individuo :: compute_fitness() { 
   double distance=0.;
   for (int i=0; i<(_ncities-1); i++) {
       distance += _distance_matrix(_x(i), _x(i + 1));
   }
   distance += _distance_matrix(_x(0), _x(_ncities - 1)); //distanza per tornare a casa
   _fitness = distance;
   return ;
}



uvec Individuo :: get_x() const{
   return _x;
}


void Individuo :: set_x(const uvec &x) {
   _x=x;
   return;
}

void Individuo :: set_distancematrix(const mat &distance_matrix){
   _distance_matrix = distance_matrix;
   return;
}






void Individuo :: pair_permutation(Random& _rnd) { 
   int i = static_cast<int>(_rnd.Rannyu(1, _ncities));
   int j = static_cast<int>(_rnd.Rannyu(1, _ncities));
    do {
        j = static_cast<int>(_rnd.Rannyu(1, _ncities));
    } while (i == j);
   this->swap(i, j);
   this->check();
   return;
}


void Individuo::shift(Random& _rnd) {
   int m = static_cast<int>(_rnd.Rannyu(2, _ncities-2)); //dimensione del blocco da shiftare
    int n = static_cast<int>(_rnd.Rannyu(1, _ncities - m )); //di quanto shiftare
    int i = static_cast<int>(_rnd.Rannyu(1, _ncities - n - m)); // indice prima città del blocco
    uvec y = _x;
    for (int k = 0; k < m; k++) {
      _x(i+n+k)=y(i+k);
    } 
    for (int k = 0; k < n; k++) {
      _x(i+k)=y(i+k+m);
    } 
   this->check();
   return;
}


void Individuo::contiguous_permutation(Random& _rnd) {
   int m = static_cast<int>(_rnd.Rannyu(2, _ncities / 2));
    int i = static_cast<int>(_rnd.Rannyu(1, _ncities - 2*m));
    int j = static_cast<int>(_rnd.Rannyu(i+m , _ncities - m + 1 ));;
    
    for (int k = 0; k < m; k++) {
      this->swap(i , j);
      i++;
      j++;
    }
   this->check();
   return;
}


void Individuo :: inversion(Random& _rnd) { //indici su ipad poi commenta se no non si capisce 

   int i = static_cast<int>(_rnd.Rannyu(1, _ncities - 3)); //indice inizio
    int l = static_cast<int>(_rnd.Rannyu(2, _ncities-i+1)); //lunghezza blocco
    int j = i+l-1; //fine
    for (int k = 0; k<l/2; k++) {
      this->swap(i, j);
      i++;
      j--;
    }
   this->check();
   return;
}


