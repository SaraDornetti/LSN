#include <iostream>
#include <cmath>
#include <algorithm>
#include <unordered_set>
#include <vector>
#include <mpi.h>
#include "system.h"

using namespace std;
using namespace arma;

// INIZIALIZZAZIONE SISTEMA
void System :: initialize(Random& _rnd) {
    
    ifstream input("input.dat"); // Start reading input.dat
    string property;
    while ( !input.eof() ){
        input >> property;
        if( property == "N_CITIES" ){
            input >> _ncities;
        }
        else if( property == "N_INDIVIDUI" ){
            input >> _nindividui;
        }
       else if( property == "CITY_CONFIGURATION" ){
            input >> _city_shape;
            if(_city_shape > 3){
                cerr << "PROBLEM: unknown city shape" << endl;
                exit(EXIT_FAILURE);
            }
        }
        else if( property == "N_GENERATIONS" ){
            input >> _ngenerations;
        }
        else if( property == "M_RATE" ){
            input >> _mutation_rate;
        }
        else if( property == "C_RATE" ){
            input >> _crossover_rate;
        }
        else if( property == "PESO" ){
            input >> _p;
        }
        else if( property == "ENDINPUT" ){
            break;
        } else cerr << "PROBLEM: unknown input" << endl;
    }
    input.close();

    _generation_index=0;

    if (_city_shape == 0) this->initialize_circle(_rnd);
    if (_city_shape == 1) this->initialize_square(_rnd);
    if (_city_shape == 2) this->initialize_Italy();

    _individuo.set_size(_nindividui);
    _individuo_figlio.set_size(_nindividui);

    for (int i=0; i<_nindividui; i++){
        _individuo(i).setncities(_ncities);
        _individuo_figlio(i).setncities(_ncities);
        _individuo(i).set_distancematrix(_distance_matrix);
        _individuo(i).initialize(_rnd);
        //_individuo_figlio(i).initialize(_rnd);
    }
    this->selection();
    return;
}


// INIZIALIZZAZIONE CERCHIO
void System :: initialize_circle(Random& _rnd) {

    vec theta_cities;
    theta_cities.set_size(_ncities); 
    ofstream coutf("configurazione_città.dat");

    if(coutf.is_open()){
        for ( int i=0; i < _ncities; i++ ) {
            theta_cities(i)=_rnd.Rannyu(0, 2*M_PI);
            coutf << i << setw(16) << theta_cities(i) << endl;
        }
        coutf.close();
    }
    else {
        cerr << "PROBLEM: Unable to open configurazione_città.dat" << endl;
    }

    mat D(_ncities, _ncities, fill::zeros);
    for ( int i=0; i<_ncities-1; i++ ) {
        for (int j=i+1; j<_ncities; j++ ) {
            D(i, j) = sqrt(2*(1-cos(theta_cities(i)-theta_cities(j))));
            D(j, i) = D(i, j);
        }
    }
    _distance_matrix = D;
    return;
 }


// INIZIALIZZAZIONE QUADRATO
void System :: initialize_square(Random& _rnd) {

    vec x_cities;
    vec y_cities;
    x_cities.set_size(_ncities); 
    y_cities.set_size(_ncities); 
    ofstream coutf("configurazione_città.dat");

    if(coutf.is_open()){
        for ( int i=0; i < _ncities; i++ ) {
            x_cities(i)=_rnd.Rannyu();
            y_cities(i)=_rnd.Rannyu();
            coutf << i << setw(16) << x_cities(i) << setw(16) << y_cities(i) << endl;
        }
        coutf.close();
    }
    else cerr << "PROBLEM: Unable to open configurazione_città.dat" << endl;

    mat D(_ncities, _ncities, fill::zeros);
    for ( int i=0; i<_ncities-1; i++ ) {
        for (int j=i+1; j<_ncities; j++ ) {
            D(i, j) = sqrt(pow(x_cities(i) - x_cities(j), 2) + pow(y_cities(i) - y_cities(j), 2));
            D(j, i) = D(i, j);
        }
    }
    _distance_matrix = D;
}


// INIZIALIZZAZIONE ITALIA
void System :: initialize_Italy() {

    vec x_cities;
    vec y_cities;
    x_cities.set_size(_ncities); 
    y_cities.set_size(_ncities); 

    ifstream cinf;
    cinf.open("cap_prov_ita.dat");
    if(cinf.is_open()){
        double x, y;
        for(int i=0; i<_ncities; i++){
            cinf >> x >> y; 
            x_cities(i)=x;
            y_cities(i)=y;
        }
    } else cerr << "PROBLEM: Unable to open INPUT file cap_prov_ita.dat"<< endl;
    cinf.close();
  
    mat D(_ncities, _ncities, fill::zeros);
    for ( int i=0; i<_ncities-1; i++ ) {
        for (int j=i+1; j<_ncities; j++ ) {
            D(i, j) = (pow(x_cities(i) - x_cities(j), 2) + pow(y_cities(i) - y_cities(j), 2));
            D(j, i) = D(i, j);
        }
    }
    _distance_matrix = D;
 }


void System :: new_generation (Random& _rnd, int rank) {

    int elite = 1;
    int n_mut = static_cast<int>(_mutation_rate*_nindividui/4);
    int n_crossover = static_cast<int>(_crossover_rate*_nindividui);
        
    int idx = elite;

    for (int i=idx; i<_nindividui-1; i += 2) {
        crossover(_rnd, i, i+1);
    }


    for (int i=0; i<n_mut; i++){
        int n = static_cast<int>(_rnd.Rannyu(elite,_nindividui));
        _individuo_figlio(n).pair_permutation(_rnd);
    }
    
    for (int i=0; i<n_mut; i++){
        int n = static_cast<int>(_rnd.Rannyu(elite,_nindividui));
        _individuo_figlio(n).shift(_rnd);
    }
    
    
    //if (_generation_index<_ngenerations*0.8){ //nella fase finale evitiamo di mescolare troppo
        for (int i=0; i<n_mut; i++){
            int n = static_cast<int>(_rnd.Rannyu(elite,_nindividui));
            _individuo_figlio(n).contiguous_permutation(_rnd);
        }
        for (int i=0; i<n_mut; i++){
            int n = static_cast<int>(_rnd.Rannyu(elite,_nindividui));
            _individuo_figlio(n).inversion(_rnd);
        }
    //}
    

    for(int i=idx; i<_nindividui; i++) {
        uvec path = _individuo_figlio(i).get_x();
        _individuo(i).set_x(path);
        _individuo(i).compute_fitness();
    }

    _generation_index ++;
    this->selection();
    this->mean_fitness(rank);

    return;
}




int System :: get_ngenerations() {
    return _ngenerations;
}


void System :: selection() {
    vector<int> indici(_nindividui);
    for (int i=0; i< _nindividui; i++){
        indici[i]=i;
    }

    // Ordina gli indici in base al fitness (fitness minore = migliore)
    sort(indici.begin(), indici.end(),
         [&](int a, int b) {
             return _individuo(a).get_fitness() < _individuo(b).get_fitness();
         });

    // Crea un nuovo field ordinato
    field<Individuo> ordinato(_nindividui);
    for (int i = 0; i < _nindividui; i++)
        ordinato(i) = _individuo(indici[i]);

    _individuo = ordinato;
    return;
}



// Algorithm for selecting the individuals:
// - p = 1: completely random selection, with no preference for the best or worst.
// - p < 1: tends to favor the poorest individuals, thus more exploration.
// - p > 1: tends to favor the best individuals, thus more exploitation. This will make best individuals reproduce more.
// The larger p is, the more likely individuals in the top positions (the best, if the population is ORIDNATED by increasing fitness) are chosen.  

//int(_nind*pow(rnd.Rannyu(), p));


void System::crossover(Random& _rnd, int ind1, int ind2) {
    int a, b;

    do {
        a = static_cast<int>(pow(_rnd.Rannyu(), _p) * (_nindividui));
        b = static_cast<int>(pow(_rnd.Rannyu(), _p) * (_nindividui));
    } while (b == a);

    uvec parent1 = _individuo(a).get_x();
    uvec parent2 = _individuo(b).get_x();

    int cut = static_cast<int>(_rnd.Rannyu(1, _ncities - 1));


    // Primo figlio
    uvec first_son(_ncities);
    first_son.fill(UINT_MAX);
    for (int i = 0; i < cut; i++) first_son(i) = parent1(i);
    unordered_set<unsigned int> used1(parent1.begin(), parent1.begin() + cut);
    int idx = cut;
    for (int i = 0; i < _ncities && idx < _ncities; i++) {
        unsigned int city = parent2(i);
        if (used1.find(city) == used1.end())
            first_son(idx++) = city;
    }

    // Secondo figlio
    uvec second_son(_ncities);
    second_son.fill(UINT_MAX);
    for (int i = 0; i < cut; i++) second_son(i) = parent2(i);
    unordered_set<unsigned int> used2(parent2.begin(), parent2.begin() + cut);
    idx = cut;
    for (int i = 0; i < _ncities && idx < _ncities; i++) {
        unsigned int city = parent1(i);
        if (used2.find(city) == used2.end())
            second_son(idx++) = city;
    }

    _individuo_figlio(ind1).set_x(first_son);
    _individuo_figlio(ind2).set_x(second_son);

    //_individuo_figlio(ind1).check();
    //_individuo_figlio(ind2).check();
    
}



void System :: write_configuration(){
    uvec best_path = _individuo(0).get_x();

    ofstream coutf("best_path.dat");
    if(coutf.is_open()){
        for ( int i=0; i < _ncities; i++ ) {
            coutf << i << setw(16) << best_path(i) << endl;
        }
        coutf.close();
    }
    else cerr << "PROBLEM: Unable to open best_path.dat" << endl;
    return;
}



void System :: finalize(){
    this->write_configuration();
    return;
}


void System :: mean_fitness(int rank) {

    double sum_fitness=0;
    for (int i=0; i<_nindividui/2; i++) {
       sum_fitness += _individuo(i).get_fitness();
    }

    double mean = sum_fitness/(_nindividui/2);

    ofstream couta;
    couta.open("mean_fitness.dat",ios::app);
    couta << rank << setw(16) << _generation_index << setw(16) << _individuo(0).get_fitness() << setw(16) << mean  << endl;
    couta.close();

    return;

}


void System :: migration(int rank, int size) {

    uvec individuo_best = _individuo(0).get_x(); //individuo da inviare
    vector<int> ind_send(_ncities);
    for (int i=0; i<_ncities; i++){
        ind_send[i] = individuo_best(i);
    }
    vector<int> ind_receive(_ncities); //vettore che conterrà indici dell'individuo da ricevere

    int send_to = (rank + 1) % size; //processo che riceve
    int recv_from = (rank - 1 + size) % size; //processo che invia

   // cout << "Rank " << rank << ": send tour to rank " << send_to << ", receive from rank " << recv_from << endl;

    MPI_Sendrecv(ind_send.data(), _ncities, MPI_INT, send_to, 0, ind_receive.data(), _ncities, MPI_INT, recv_from, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

   // cout << "Rank " << rank << ": ricevuto tour: ";
   // for (int c : ind_receive) cout << c << " ";
   // cout << endl;
  
    uvec new_individuo_best(_ncities); 
    for (int i=0; i<_ncities; i++){
        new_individuo_best(i) = ind_receive[i]; 
    }

    _individuo(0).set_x(new_individuo_best);
    _individuo(0).compute_fitness();

}


Individuo System :: get_individuo(int idx) {
    return _individuo(idx);
}