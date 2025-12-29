/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <cmath>
#include <cstdlib>
#include <string>
#include "system.h"

using namespace std;
using namespace arma;

void System :: step(){ // Perform a simulation step
  if(_sim_type == 0) this->Verlet();  // Perform a MD step
  else for(int i=0; i<_npart; i++) this->move(int(_rnd.Rannyu()*_npart)); // Perform a MC step on a randomly choosen particle
  _nattempts += _npart; //update number of attempts performed on the system
  return;
}

void System :: Verlet(){
  double xnew, ynew, znew;
  for(int i=0; i<_npart; i++){ //Force acting on particle i
    _fx(i) = this->Force(i,0); // componente lungo x della forza agente su particelle i-esima
    _fy(i) = this->Force(i,1); // componente lungo y della forza agente su particelle i-esima
    _fz(i) = this->Force(i,2); // componente lungo z della forza agente su particelle i-esima
  }
  for(int i=0; i<_npart; i++){ //Verlet integration scheme
    xnew = this->pbc( 2.0 * _particle(i).getposition(0,true) - _particle(i).getposition(0,false) + _fx(i) * pow(_delta,2), 0); // calcolo la nuova x usando lo schema Verlet, prestando attenzione alle pbc
    ynew = this->pbc( 2.0 * _particle(i).getposition(1,true) - _particle(i).getposition(1,false) + _fy(i) * pow(_delta,2), 1); // calcolo la nuova y usando lo schema Verlet, prestando attenzione alle pbc 
    znew = this->pbc( 2.0 * _particle(i).getposition(2,true) - _particle(i).getposition(2,false) + _fz(i) * pow(_delta,2), 2); // calcolo la nuova z usando lo schema Verlet, prestando attenzione alle pbc
    //getposition(0, true) prende la posizione attuale lungo x, false quella vecchia lungo x (0->x, 1->y, 2->z)
    _particle(i).setvelocity(0, this->pbc(xnew - _particle(i).getposition(0,false), 0)/(2.0 * _delta));
    _particle(i).setvelocity(1, this->pbc(ynew - _particle(i).getposition(1,false), 1)/(2.0 * _delta));
    _particle(i).setvelocity(2, this->pbc(znew - _particle(i).getposition(2,false), 2)/(2.0 * _delta));
    _particle(i).acceptmove(); // xold = xnew
    //penso che intenda che copia la posizione attuale in quella vecchia, nella riga commentata sopra xnew non sono quelle appena calcolate
    _particle(i).setposition(0, xnew);
    _particle(i).setposition(1, ynew);
    _particle(i).setposition(2, znew);
  }
  _naccepted += _npart; //in MD non si accetta o rifiuta come nel metropolis, quindi il numero di mosse accettate è uguale al numero di mosse fatte
  return;
}

void System :: ReverseVerlet(){
  if (!_reverse) return;

  for (int i = 0; i < _npart; i++) {
    
    //SALVO LE POSIZIONI LETTE DA FILE
    // posizioni correnti salvate in conf.xyz
    double current_x = _particle(i).getposition(0, true);   
    double current_y = _particle(i).getposition(1, true);
    double current_z = _particle(i).getposition(2, true);

    // posizioni allo step precedente salvate in conf-1.xyz
    double old_x = _particle(i).getposition(0, false);  
    double old_y = _particle(i).getposition(1, false);
    double old_z = _particle(i).getposition(2, false);

    //INVERTO POSIZIONI CORRENTI E POSIZIONI VECCHIE PER IMPLEMENTARE REVERSE
    // salvo le posizioni correnti e accetto la mossa -> in questo modo le posizioni correnti diventano le posizioni vecchie
    _particle(i).setposition(0, current_x);
    _particle(i).setposition(1, current_y);
    _particle(i).setposition(2, current_z);
    _particle(i).acceptmove();

    // salvo le posizioni vecchie come correnti
    _particle(i).setposition(0, old_x);
    _particle(i).setposition(1, old_y);
    _particle(i).setposition(2, old_z);

    //CALCOLO LE VELOCITA' INIZIALI CON IL SEGNO CORRETTO PER TORNARE INDIETRO
    _particle(i).setvelocity(0, this->pbc(old_x - current_x, 0) / (2.0 * _delta));
    _particle(i).setvelocity(1, this->pbc(old_y - current_y, 1) / (2.0 * _delta));
    _particle(i).setvelocity(2, this->pbc(old_z - current_z, 2) / (2.0 * _delta));
  }

  _reverse = false;
}


double System :: Force(int i, int dim){ //calcola la componente lungo dim della forza (Lennard-Jones) agente sulla particella i-esima, usando un cutoff sferico _r_cut
  double f=0.0, dr;
  vec distance;
  distance.resize(_ndim); //in questo modo il vettore distanza ha lo stesso numero di componenti delle dimensioni del sistema
  for (int j=0; j<_npart; j++){
    if(i != j){ //calcola la distanza tra la particella i-esima e tutte le altre, quindi esclude la stessa
      distance(0) = this->pbc( _particle(i).getposition(0,true) - _particle(j).getposition(0,true), 0); 
      distance(1) = this->pbc( _particle(i).getposition(1,true) - _particle(j).getposition(1,true), 1);
      distance(2) = this->pbc( _particle(i).getposition(2,true) - _particle(j).getposition(2,true), 2);
      //posizioni attuali, lungo x, y e z.
      dr = sqrt( dot(distance,distance) ); //dot sarà una funzione di armadillo che implementa il prodotto scalare tra due vettori
      if(dr < _r_cut){ //approssimazione per cui consideriamo forze soltanto entro un certo range di distanza
        f += distance(dim) * (48.0/pow(dr,14) - 24.0/pow(dr,8)); //non lo azzero perchè lo faccio solo per la particella i-esima, quindi se richiamo la funzione riazzera f
      //forza lungo una componente (dal potenziale di Lennard-Jones)
      }
    }
  }
  return f;
}


// a seconda del tipo di simulazione (_sim_type) cambia il comportamento della funzione.
void System :: move(int i){ // Propose a MC move for particle i !!!!
  if(_sim_type == 3){ //Gibbs sampler for Ising
    // TO BE FIXED IN EXERCISE 6
  } else {           // M(RT)^2
    if(_sim_type == 1){       // LJ system
      vec shift(_ndim);       // Will store the proposed translation //si costruisce un vettore della stessa dimensione del sistema
      for(int j=0; j<_ndim; j++){
        shift(j) = _rnd.Rannyu(-1.0,1.0) * _delta; // uniform distribution in [-_delta;_delta)
      }
      _particle(i).translate(shift, _side);  //Call the function Particle::translate //sposto la particella i-esima
      if(this->metro(i)){ //Metropolis acceptance evaluation
        _particle(i).acceptmove();
        _naccepted++;
      } else _particle(i).moveback(); //If translation is rejected, restore the old configuration
    } else {                  // Ising 1D
      if(this->metro(i)){     //Metropolis acceptance evaluation for a spin flip involving spin i
        _particle(i).flip();  //If accepted, the spin i is flipped
        _naccepted++;
      }
    }
  }
  return;
}

bool System :: metro(int i){ // Metropolis algorithm
  bool decision = false;
  double delta_E, acceptance;
  if(_sim_type == 1) delta_E = this->Boltzmann(i,true) - this->Boltzmann(i,false);
  else delta_E = 2.0 * _particle(i).getspin() * 
                 ( _J * (_particle(this->pbc(i-1)).getspin() + _particle(this->pbc(i+1)).getspin() ) + _H );
  acceptance = exp(-_beta*delta_E);
  if(_rnd.Rannyu() < acceptance ) decision = true; //Metropolis acceptance step
  return decision;
}

//Boltzmann serve per calcolare il peso di Boltzmann per lo step metropolis
double System :: Boltzmann(int i, bool xnew){ //Calcola energia potenziale di una particella i con la configurazione scelta (nuova o vecchia)
  double energy_i=0.0;
  double dx, dy, dz, dr;
  for (int j=0; j<_npart; j++){
    if(j != i){ //considero tutte le particelle
      //per la particella i-esima inserisco nella funzione se voglio l'energia vecchia o nuova (bool xnew)
      //per le particelle j-esime prende la configurazione 1
      dx = this->pbc(_particle(i).getposition(0,xnew) - _particle(j).getposition(0,1), 0);
      dy = this->pbc(_particle(i).getposition(1,xnew) - _particle(j).getposition(1,1), 1);
      dz = this->pbc(_particle(i).getposition(2,xnew) - _particle(j).getposition(2,1), 2);
      dr = dx*dx + dy*dy + dz*dz; //calcolo il quadrato della distanza
      dr = sqrt(dr); //la distanza
      if(dr < _r_cut){ //considero solo particelle entro una sfera di raggio _r_cut
        energy_i += 1.0/pow(dr,12) - 1.0/pow(dr,6); //calcolo il potenziale di LJ (manca fattore 4 qui)
      }
    }
  }
  return 4.0 * energy_i; //aggiungo il fattore 4
}

double System :: pbc(double position, int i){ // Enforce periodic boundary conditions
//position: la coordinata (ad esempio x, y o z) della particella. Può assumere qualsiasi valore reale (anche al di fuori della scatola simulata).
//i: l'indice della direzione (0 = x, 1 = y, 2 = z, ecc.), serve per sapere la lunghezza della scatola in quella direzione.
  return position - _side(i) * rint(position / _side(i));
  //Divide la coordinata per la lunghezza della scatola -> ottiene quanti "box" di lato _side(i) la particella ha superato.
  //rint() arrotonda al più vicino intero (round to nearest).
  //Sottraendo questo multiplo della lunghezza della scatola, riporti la particella dentro la scatola centrale.
}

int System :: pbc(int i){ // Enforce periodic boundary conditions for spins
  //Qui i è l'indice di un sito/spin della catena.
  //_npart = numCon le PBC la catena è chiusa ad anello: lo spin subito dopo l'ultimo è il primo, e quello prima del primo è l'ultimo.
  if(i >= _npart) i = i - _npart;
  else if(i < 0)  i = i + _npart;
  return i;
} 

void System :: initialize(){ // Initialize the System object according to the content of the input files in the ../INPUT/ directory

  int p1, p2; // Read from ../INPUT/Primes a pair of numbers to be used to initialize the RNG
  ifstream Primes("../INPUT/Primes");
  Primes >> p1 >> p2 ; //Legge due numeri primi da Primes -> servono a costruire la sequenza casuale.
  Primes.close();
  int seed[4]; // Read the seed of the RNG
  ifstream Seed("../INPUT/seed.in");
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3]; //Legge 4 numeri da seed.in -> seed iniziale.
  _rnd.SetRandom(seed,p1,p2); //Insieme inizializzano il generatore di numeri casuali _rnd.

  ofstream couta("../OUTPUT/acceptance.dat"); // Set the heading line in file ../OUTPUT/acceptance.dat
  couta << "#   N_BLOCK:  ACCEPTANCE:" << endl;
  couta.close(); //Crea (o resetta) acceptance.dat, con intestazione per memorizzare le accettazioni delle mosse Monte Carlo.

  //Apre input.dat per leggere i parametri della simulazione.
  // Apre output.dat per scrivere le informazioni iniziali.
  ifstream input("../INPUT/input.dat"); // Start reading ../INPUT/input.dat
  ofstream coutf;
  coutf.open("../OUTPUT/output.dat"); 
  string property;
  double delta;
  while ( !input.eof() ){
    input >> property; //Legge parola chiave (property) -> decide quale variabile del sistema settare.
    if( property == "SIMULATION_TYPE" ){
      input >> _sim_type;
      if(_sim_type > 1){
        input >> _J;
        input >> _H;
      }
      if(_sim_type > 3){
        cerr << "PROBLEM: unknown simulation type" << endl;
        exit(EXIT_FAILURE);
      }
      if(_sim_type == 0)      coutf << "LJ MOLECULAR DYNAMICS (NVE) SIMULATION"  << endl;
      else if(_sim_type == 1) coutf << "LJ MONTE CARLO (NVT) SIMULATION"         << endl;
      else if(_sim_type == 2) coutf << "ISING 1D MONTE CARLO (MRT^2) SIMULATION" << endl;
      else if(_sim_type == 3) coutf << "ISING 1D MONTE CARLO (GIBBS) SIMULATION" << endl;
    } else if( property == "RESTART" ){
      input >> _restart;
    } else if( property == "TEMP" ){
      input >> _temp;
      _beta = 1.0/_temp;
      coutf << "TEMPERATURE= " << _temp << endl; //salva in _temp, calcola beta = 1/T.
    } else if( property == "NPART" ){
      input >> _npart;
      //Alloca array di forze _fx, _fy, _fz.
      _fx.resize(_npart);
      _fy.resize(_npart);
      _fz.resize(_npart);
      _particle.set_size(_npart); //Alloca array di particelle _particle.
      for(int i=0; i<_npart; i++){ 
        _particle(i).initialize(); //Inizializza ogni particella.
        if(_rnd.Rannyu() > 0.5) _particle(i).flip(); // to randomize the spin configuration
        //Per l'Ising, randomizza spin con probabilità 50%.
      }
      coutf << "NPART= " << _npart << endl;
    } else if( property == "RHO" ){
      input >> _rho;
      _volume = _npart/_rho; //Calcola il volume totale
      _side.resize(_ndim);
      _halfside.resize(_ndim);
      double side = pow(_volume, 1.0/3.0); //Lato del box
      for(int i=0; i<_ndim; i++) _side(i) = side; //box cubico
      _halfside=0.5*_side;
      coutf << "SIDE= ";
      for(int i=0; i<_ndim; i++){
        coutf << setw(12) << _side[i]; //setw(12) -> significa set width 12, cioè ogni numero verrà stampato occupando 12 spazi.
      }
      coutf << endl;
    } else if( property == "R_CUT" ){
      input >> _r_cut;
      coutf << "R_CUT= " << _r_cut << endl;
    } else if( property == "DELTA" ){
      input >> delta;
      coutf << "DELTA= " << delta << endl;
      _delta = delta;
    } else if( property == "NBLOCKS" ){
      input >> _nblocks;
      coutf << "NBLOCKS= " << _nblocks << endl;
    } else if( property == "NSTEPS" ){
      input >> _nsteps;
      coutf << "NSTEPS= " << _nsteps << endl;
    } else if( property == "ENDINPUT" ){
      coutf << "Reading input completed!" << endl;
      break;
    } else cerr << "PROBLEM: unknown input" << endl;
  }
  input.close();
  this->read_configuration(); //Legge configurazioni iniziali da file.
  if(_sim_type==0) this->initialize_velocities(); //Se è dinamica molecolare (tipo 0), inizializza anche le velocità delle particelle.
  this->ReverseVerlet();
  coutf << "System initialized!" << endl;
  coutf.close();
  return;
}


//serve per inizializzare le velocità delle molecole nel caso di dinamica molecolare (_sim_type==0)
void System :: initialize_velocities(){
  double xold, yold, zold;
  if(_restart){ //considero il caso in cui _restart!=0, quindi ricomincio da simulazione passata
    ifstream cinf;
    cinf.open("../INPUT/CONFIG/conf-1.xyz"); //Apre il file conf-1.xyz, che contiene una configurazione salvata (con posizioni old).
    if(cinf.is_open()){ 
      string comment;
      string particle;
      int ncoord;
      cinf >> ncoord;
      if (ncoord != _npart){ 
        //Controlla che il numero di particelle letto dal file (ncoord) coincida con _npart (il valore da input.dat).
        cerr << "PROBLEM: conflicting number of coordinates in input.dat & config.xyz not match!" << endl;
        exit(EXIT_FAILURE);
      }
      cinf >> comment;
      for(int i=0; i<_npart; i++){
        cinf >> particle >> xold >> yold >> zold; // units of coordinates in conf.xyz is _side
        //salva le coordinate come posizioni vecchie!!
        _particle(i).setpositold(0, this->pbc(_side(0)*xold, 0));
        _particle(i).setpositold(1, this->pbc(_side(1)*yold, 1));
        _particle(i).setpositold(2, this->pbc(_side(2)*zold, 2));
      }
    } else cerr << "PROBLEM: Unable to open INPUT file conf-1.xyz"<< endl;
    cinf.close();
  } else { // TO BE TEMPORARILY CHANGED IN EXERCISE 4
    //Se non è un restart, bisogna generare velocità nuove 
    vec vx(_npart), vy(_npart), vz(_npart); //vettori di dimensione pari al numero di particelle dove scrivo le velocità lungo x y e z di ciascuna particella
    vec sumv(_ndim); //creo un vettore di componenti pari alla dimensione del sistema (e.g. 3 componenti)
    sumv.zeros(); //inizializzo a zero

    //da qui inizia la parte nuova scritta da me, seguendo la consegna del 04.2
    vx.zeros();
    vy.zeros();
    vz.zeros();

    double v_star = sqrt(2*_temp); 
    for (int i=0; i<_npart; i++) {
      int direzione = static_cast<int>(_rnd.Rannyu(0,3));
      int segno;

      if (_rnd.Rannyu(-1,1)>0) segno = 1;
      else segno = -1;

      if (direzione == 0) vx(i) = segno*v_star;
      else if (direzione == 1) vy(i) = segno*v_star;
      else vz(i) = segno*v_star;

      sumv(0) += vx(i);
      sumv(1) += vy(i);
      sumv(2) += vz(i);

    }
    //qui finisce la parte aggiunta da me

    /*for (int i=0; i<_npart; i++){
      //genera velocità distribuite secondo una gaussiana con media 0 e dev st pari a radice della temperatura
      vx(i) = _rnd.Gauss(0.,sqrt(_temp)); 
      vy(i) = _rnd.Gauss(0.,sqrt(_temp));
      vz(i) = _rnd.Gauss(0.,sqrt(_temp));
      //calcola la somma delle velocità che serve per rimuovere il moto del centro di massa
      sumv(0) += vx(i);
      sumv(1) += vy(i);
      sumv(2) += vz(i);
    }*/

    for (int idim=0; idim<_ndim; idim++) sumv(idim) = sumv(idim)/double(_npart); //calcola la velocità media per ogni direzione
    double sumv2 = 0.0, scalef; //sumv2 è la somma di tutte le componenti di velocità al quadrato (di tutte le particelle) diviso per il numero di particelle 
    //scalef è un fattore di scala
    //calcolo le velocità togliendo quella del centro di massa
    for (int i=0; i<_npart; i++){
      vx(i) = vx(i) - sumv(0);
      vy(i) = vy(i) - sumv(1);
      vz(i) = vz(i) - sumv(2);
      //sottraendo la velocità media, ora la somma di tutte le velocità è zero
      sumv2 += vx(i) * vx(i) + vy(i) * vy(i) + vz(i) * vz(i);
    }
    sumv2 /= double(_npart); //divido per il numero di particelle -> energia cinetica media
    //riscalatura alla temperatura desiderata -> così ho energia cinetica corrispondente a T
    //in 3D energia cinetica media per particella deve essere (3/2)T
    scalef = sqrt(3.0 * _temp / sumv2);   // velocity scale factor 
    for (int i=0; i<_npart; i++){
      _particle(i).setvelocity(0, vx(i)*scalef);
      _particle(i).setvelocity(1, vy(i)*scalef);
      _particle(i).setvelocity(2, vz(i)*scalef);
    }
    for (int i=0; i<_npart; i++){ //Calcola le posizioni old usando un passo indietro di Verlet:
      xold = this->pbc( _particle(i).getposition(0,true) - _particle(i).getvelocity(0)*_delta, 0);
      yold = this->pbc( _particle(i).getposition(1,true) - _particle(i).getvelocity(1)*_delta, 1);
      zold = this->pbc( _particle(i).getposition(2,true) - _particle(i).getvelocity(2)*_delta, 2);
      _particle(i).setpositold(0, xold);
      _particle(i).setpositold(1, yold);
      _particle(i).setpositold(2, zold);
    }
  }
  return;
}


//Legge da properties.dat quali grandezze misurare (energia, pressione, g(r), ecc.)
//Crea file di output corrispondenti.
//Alloca e azzera i vettori per accumulare medie e varianze a blocchi.
void System :: initialize_properties(){ // Initialize data members used for measurement of properties

  string property; //property conterrà le parole chiave lette da properties.dat.
  int index_property = 0; //index_property è un contatore per assegnare un "indice" a ogni proprietà monitorata.
  _nprop = 0; //_nprop = numero totale di proprietà che verranno misurate (serve per dimensionare i vettori statistici).

  _measure_penergy  = false; //Defining which properties will be measured
  _measure_kenergy  = false;
  _measure_tenergy  = false;
  _measure_pressure = false;
  _measure_gofr     = false;
  _measure_magnet   = false;
  _measure_cv       = false;
  _measure_chi      = false;
  _measure_pofv     = false;
  //Tutti i flag booleani sono messi a false -> di default non si misura niente.

  ifstream input("../INPUT/properties.dat"); //Apre properties.dat
  if (input.is_open()){
    // Legge una parola alla volta -> in base al valore entra in uno specifico blocco if.
    while ( !input.eof() ){
      input >> property;
      if( property == "POTENTIAL_ENERGY" ){
        ofstream coutp("../OUTPUT/potential_energy.dat"); //per ogni file di output crea un file di output relativo 
        coutp << "#     BLOCK:  ACTUAL_PE:     PE_AVE:      ERROR:" << endl; //con intestazione
        coutp.close();
        _nprop++; //aggiorno indice proprietà
        _index_penergy = index_property; //metto  _index_penergy a 0 -> se c'è flag questa sarà la prima misura (con indice zero)
        _measure_penergy = true;
        index_property++; //aggiorno index_property (N.B. se non c'è questa flag non si aggiorna e la prima misura sarà un'altra)
        _vtail = 0.0; // TO BE FIXED IN EXERCISE 7
      } else if( property == "KINETIC_ENERGY" ){
        ofstream coutk("../OUTPUT/kinetic_energy.dat");
        coutk << "#     BLOCK:   ACTUAL_KE:    KE_AVE:      ERROR:" << endl;
        coutk.close();
        _nprop++;
        _measure_kenergy = true;
        _index_kenergy = index_property;
        index_property++;
      } else if( property == "TOTAL_ENERGY" ){
        ofstream coutt("../OUTPUT/total_energy.dat");
        coutt << "#     BLOCK:   ACTUAL_TE:    TE_AVE:      ERROR:" << endl;
        coutt.close();
        _nprop++;
        _measure_tenergy = true;
        _index_tenergy = index_property;
        index_property++;
      } else if( property == "TEMPERATURE" ){
        ofstream coutte("../OUTPUT/temperature.dat");
        coutte << "#     BLOCK:   ACTUAL_T:     T_AVE:       ERROR:" << endl;
        coutte.close();
        _nprop++;
        _measure_temp = true;
        _index_temp = index_property;
        index_property++;
      } else if( property == "PRESSURE" ){
        ofstream coutpr("../OUTPUT/pressure.dat");
        coutpr << "#     BLOCK:   ACTUAL_P:     P_AVE:       ERROR:" << endl;
        coutpr.close();
        _nprop++;
        _measure_pressure = true;
        _index_pressure = index_property;
        index_property++;
        _ptail = 0.0; // TO BE FIXED IN EXERCISE 7
      } else if( property == "GOFR" ){
        ofstream coutgr("../OUTPUT/gofr.dat");
        coutgr << "# DISTANCE:     AVE_GOFR:        ERROR:" << endl;
        coutgr.close();
        input>>_n_bins; //legge il numero di bin da file
        _nprop+=_n_bins; //Incrementa _nprop di _n_bins (perché ogni bin è una “proprietà” separata)
        _bin_size = (_halfside.min() )/(double)_n_bins;
        _measure_gofr = true;
        _index_gofr = index_property;
        index_property+= _n_bins; 
      } else if( property == "POFV" ){
        if(_sim_type > 0){
          cerr << "PROBLEM: DOES NOT MAKE SENSE COMPUTING POFV FOR MC" << endl;
          exit(EXIT_FAILURE);
        }
        ofstream coutpv("../OUTPUT/pofv.dat");
        coutpv << "# VELOCITY:     ACTUAL_V:     AVE_POFV:        ERROR:" << endl;
        coutpv.close();
        input>>_n_bins_v; //legge da file
        _nprop += _n_bins_v;
        //Calcola la larghezza di bin per le velocità.
        double vmax = 4.0 * sqrt(_temp); // cutoff ragionevole (circa 3 sigma della Maxwell-Boltzmann)
        _bin_size_v = vmax / double(_n_bins_v);
        //_bin_size_v = 4.0*_temp/(double)_n_bins_v; // TO BE FIXED IN EXERCISE 4
        _measure_pofv = true;
        _index_pofv = index_property;
        index_property += _n_bins_v;
      } else if( property == "MAGNETIZATION" ){
        ofstream coutpr("../OUTPUT/magnetization.dat");
        coutpr << "#     BLOCK:   ACTUAL_M:     M_AVE:       ERROR:" << endl;
        coutpr.close();
        _nprop++;
        _measure_magnet = true;
        _index_magnet = index_property;
        index_property++;
      } else if( property == "SPECIFIC_HEAT" ){
        ofstream coutpr("../OUTPUT/specific_heat.dat");
        coutpr << "#     BLOCK:   ACTUAL_CV:    CV_AVE:      ERROR:" << endl;
        coutpr.close();
        _nprop++;
        _measure_cv = true;
        _index_cv = index_property;
        index_property++;
      } else if( property == "SUSCEPTIBILITY" ){
        ofstream coutpr("../OUTPUT/susceptibility.dat");
        coutpr << "#     BLOCK:   ACTUAL_X:     X_AVE:       ERROR:" << endl;
        coutpr.close();
        _nprop++;
        _measure_chi = true;
        _index_chi = index_property;
        index_property++;
      } else if( property == "ENDPROPERTIES" ){
        ofstream coutf;
        coutf.open("../OUTPUT/output.dat",ios::app);
        coutf << "Reading properties completed!" << endl;
        coutf.close();
        break;
      } else cerr << "PROBLEM: unknown property" << endl;
    }
    input.close();
  } else cerr << "PROBLEM: Unable to open properties.dat" << endl;

  // according to the number of properties, resize the vectors _measurement,_average,_block_av,_global_av,_global_av2
  _measurement.resize(_nprop); //valori istantanei
  _average.resize(_nprop); //media sui blocchi
  _block_av.resize(_nprop); //media corrente nel blocco
  _global_av.resize(_nprop); //media su tutti i blocchi
  _global_av2.resize(_nprop); //quadrati delle medie
  //mette a zero i vettori di medie
  _average.zeros();
  _global_av.zeros();
  _global_av2.zeros();
  _nattempts = 0;
  _naccepted = 0;
  return;
}


//Scrive configurazione finale (in ../OUTPUT/CONFIG/)
//Salva lo stato del generatore random
//Scrive messaggio di completamento.
void System :: finalize(){
  this->write_configuration();
  _rnd.SaveSeed();
  ofstream coutf;
  coutf.open("../OUTPUT/output.dat",ios::app);
  coutf << "Simulation completed!" << endl;
  coutf.close();
  return;
}

// Write final configurations as .xyz files in the directory ../OUTPUT/CONFIG/
void System :: write_configuration(){
  ofstream coutf;
  if(_sim_type < 2){
    coutf.open("../OUTPUT/CONFIG/config.xyz");
    if(coutf.is_open()){
      coutf << _npart << endl;
      coutf << "#Comment!" << endl;
      for(int i=0; i<_npart; i++){
        coutf << "LJ" << "  " 
              << setprecision(17) << _particle(i).getposition(0,true)/_side(0) << "   " // x
              << setprecision(17) << _particle(i).getposition(1,true)/_side(1) << "   " // y
              << setprecision(17) << _particle(i).getposition(2,true)/_side(2) << endl; // z
      }
    } else cerr << "PROBLEM: Unable to open config.xyz" << endl;
    coutf.close();
    coutf.open("../OUTPUT/CONFIG/conf-1.xyz");
    if(coutf.is_open()){
      coutf << _npart << endl;
      coutf << "#Comment!" << endl;
      for(int i=0; i<_npart; i++){
        coutf << "LJ" << "  " //"LJ" = etichetta per particella Lennard-Jones
              << setprecision(17) << _particle(i).getposition(0,false)/_side(0) << "   " // x
              << setprecision(17) << _particle(i).getposition(1,false)/_side(1) << "   " // y
              << setprecision(17) << _particle(i).getposition(2,false)/_side(2) << endl; // z
      }
    } else cerr << "PROBLEM: Unable to open conf-1.xyz" << endl;
    coutf.close();
  } else {
    coutf.open("../OUTPUT/CONFIG/config.spin");
    for(int i=0; i<_npart; i++) coutf << _particle(i).getspin() << " ";
    coutf.close();
  }
  return;
}

// Write configuration nconf as a .xyz file in directory ../OUTPUT/CONFIG/
void System :: write_XYZ(int nconf){
  ofstream coutf;
  coutf.open("../OUTPUT/CONFIG/config_" + to_string(nconf) + ".xyz");
  if(coutf.is_open()){
    coutf << _npart << endl;
    coutf << "#Comment!" << endl;
    for(int i=0; i<_npart; i++){
      coutf << "LJ" << "  " 
            << setw(16) << _particle(i).getposition(0,true)          // x
            << setw(16) << _particle(i).getposition(1,true)          // y
            << setw(16) << _particle(i).getposition(2,true) << endl; // z
    }
  } else cerr << "PROBLEM: Unable to open config.xyz" << endl;
  coutf.close();
  return;
}

// Read configuration from a .xyz file in directory ../OUTPUT/CONFIG/
void System :: read_configuration(){
  ifstream cinf;
  cinf.open("../INPUT/CONFIG/config.xyz"); //Apre il file config.xyz, che contiene le coordinate iniziali delle particelle.
  if(cinf.is_open()){
    string comment;
    string particle;
    double x, y, z;
    int ncoord; //ncoord: numero di particelle letto dal file.
    cinf >> ncoord;
    if (ncoord != _npart){ //Viene controllato che sia uguale a _npart (definito in input.dat)
      cerr << "PROBLEM: conflicting number of coordinates in input.dat & config.xyz not match!" << endl;
      exit(EXIT_FAILURE);
    }
    cinf >> comment;
    for(int i=0; i<_npart; i++){
      cinf >> particle >> x >> y >> z; // units of coordinates in conf.xyz is _side
      _particle(i).setposition(0, this->pbc(_side(0)*x, 0)); //qui _halfside all'inzio.
      _particle(i).setposition(1, this->pbc(_side(1)*y, 1));
      _particle(i).setposition(2, this->pbc(_side(2)*z, 2));
      _particle(i).acceptmove(); // _x_old = _x_new
      //significa che la posizione "nuova" diventa anche la "vecchia" -> quindi la configurazione iniziale è accettata come stato di partenza.
    }
  } else cerr << "PROBLEM: Unable to open INPUT file config.xyz"<< endl;
  cinf.close();
  if(_restart and _sim_type > 1){
    int spin;
    cinf.open("../INPUT/CONFIG/config.spin");
    for(int i=0; i<_npart; i++){
      cinf >> spin;
      _particle(i).setspin(spin);
    }
    cinf.close();
  }
  return;
}

void System :: block_reset(int blk){ // Reset block accumulators to zero
  ofstream coutf;
  if(blk>0){
    coutf.open("../OUTPUT/output.dat",ios::app); //L’apertura è fatta con ios::app -> modalità append, cioè non sovrascrive il file ma aggiunge in fondo.
    coutf << "Block completed: " << blk << endl; //dove N è il numero del blocco appena concluso.
    coutf.close();
  }
  _block_av.zeros(); //_block_av è un vettore che tiene i valori medi delle proprietà (energia, temperatura, ecc.) all’interno del blocco.
  return;
}

void System :: measure(){ // Measure properties
  _measurement.zeros(); //Azzera il vettore _measurement, che conterrà i valori misurati in questo step.
  // POTENTIAL ENERGY, VIRIAL, GOFR ///////////////////////////////////////////
  int bin;
  vec distance;
  distance.resize(_ndim);
  double penergy_temp=0.0, dr; // temporary accumulator for potential energy
  double kenergy_temp=0.0; // temporary accumulator for kinetic energy
  double tenergy_temp=0.0;
  double magnetization=0.0;
  double virial=0.0;
  if (_measure_penergy or _measure_pressure or _measure_gofr) {
    for (int i=0; i<_npart-1; i++){
      for (int j=i+1; j<_npart; j++){
        distance(0) = this->pbc( _particle(i).getposition(0,true) - _particle(j).getposition(0,true), 0);
        distance(1) = this->pbc( _particle(i).getposition(1,true) - _particle(j).getposition(1,true), 1);
        distance(2) = this->pbc( _particle(i).getposition(2,true) - _particle(j).getposition(2,true), 2);
        dr = sqrt( dot(distance,distance) );
        // GOFR ... TO BE FIXED IN EXERCISE 7
        if(dr < _r_cut){
          if(_measure_penergy)  penergy_temp += 1.0/pow(dr,12) - 1.0/pow(dr,6); // POTENTIAL ENERGY
          if(_measure_pressure) virial       += 1.0/pow(dr,12) - 0.5/pow(dr,6); // PRESSURE
        }
      }
    }
  }

  // POFV ... TO BE FIXED IN EXERCISE 4
  if (_measure_pofv){ // _measure_pofv flag for measuring the velocity modulus distribution
    for (int i=0; i<_npart; i++) { //devo capire come fare la media sui blocchi
      double v_mod = sqrt(dot( _particle(i).getvelocity() , _particle(i).getvelocity() )); //calcolo il modulo delle velocità per ciascuna particella
      int bin = int(v_mod/_bin_size_v); //assegno alla velocità un bin
      if(bin<_n_bins_v){ //controllo
        _measurement(_index_pofv+bin) += 1.0/(double)_npart; //normalizzato sul numero di particelle
      }
    }
  }

  // POTENTIAL ENERGY //////////////////////////////////////////////////////////
  if (_measure_penergy){
    penergy_temp = _vtail + 4.0 * penergy_temp / double(_npart); //Aggiunge la correzione del tail (_vtail), che compensa il troncamento del potenziale.
    _measurement(_index_penergy) = penergy_temp;
  }
  // KINETIC ENERGY ////////////////////////////////////////////////////////////
  if (_measure_kenergy){
    for (int i=0; i<_npart; i++) kenergy_temp += 0.5 * dot( _particle(i).getvelocity() , _particle(i).getvelocity() ); 
    kenergy_temp /= double(_npart);
    _measurement(_index_kenergy) = kenergy_temp;
  }
  // TOTAL ENERGY (kinetic+potential) //////////////////////////////////////////
  if (_measure_tenergy){
    if (_sim_type < 2) _measurement(_index_tenergy) = kenergy_temp + penergy_temp;
    else { 
      double s_i, s_j;
      for (int i=0; i<_npart; i++){
        s_i = double(_particle(i).getspin());
        s_j = double(_particle(this->pbc(i+1)).getspin());
        tenergy_temp += - _J * s_i * s_j - 0.5 * _H * (s_i + s_j);
      }
      tenergy_temp /= double(_npart);
      _measurement(_index_tenergy) = tenergy_temp;
    }
  }
  // TEMPERATURE ///////////////////////////////////////////////////////////////
  if (_measure_temp and _measure_kenergy) _measurement(_index_temp) = (2.0/3.0) * kenergy_temp;
  // PRESSURE //////////////////////////////////////////////////////////////////
  if (_measure_pressure) _measurement[_index_pressure] = _rho * (2.0/3.0) * kenergy_temp + (_ptail*_npart + 48.0*virial/3.0)/_volume; //include correzione del tail _ptail.
  // MAGNETIZATION /////////////////////////////////////////////////////////////
// TO BE FIXED IN EXERCISE 6
  // SPECIFIC HEAT /////////////////////////////////////////////////////////////
// TO BE FIXED IN EXERCISE 6
  // SUSCEPTIBILITY ////////////////////////////////////////////////////////////
// TO BE FIXED IN EXERCISE 6

//Aggiorna il vettore _block_av (somme all'interno del blocco).
  _block_av += _measurement; //Update block accumulators

  return;
}

void System :: averages(int blk){

  ofstream coutf;
  double average, sum_average, sum_ave2;

  _average     = _block_av / double(_nsteps); 
  //_block_av = accumulatore con la somma delle osservabili nel blocco corrente.
  //Dividendo per _nsteps ottieni la media del blocco (_average).
  //_nsteps=// Number of simulation steps in each block
  _global_av  += _average; //_global_av = somma cumulativa delle medie dei blocchi (serve per calcolare la media globale).
  _global_av2 += _average % _average; // % -> element-wise multiplication
  //_global_av2 = somma cumulativa dei quadrati delle medie dei blocchi (serve per calcolare la varianza).
  //Il simbolo % non è l'operatore modulo qui: è sovraccaricato per fare la moltiplicazione elemento per elemento 
  //tra vettori (element-wise). Quindi _average % _average = il vettore dei quadrati di _average.

  // POTENTIAL ENERGY //////////////////////////////////////////////////////////
  if (_measure_penergy){
    coutf.open("../OUTPUT/potential_energy.dat",ios::app);
    average  = _average(_index_penergy); // media del blocco
    sum_average = _global_av(_index_penergy); // somma globale
    sum_ave2 = _global_av2(_index_penergy); // somma quadrati
    coutf << setw(12) << blk 
          << setw(12) << average // valore del blocco
          << setw(12) << sum_average/double(blk) // media cumulativa
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  // KINETIC ENERGY ////////////////////////////////////////////////////////////
  if (_measure_kenergy){
    coutf.open("../OUTPUT/kinetic_energy.dat",ios::app);
    average  = _average(_index_kenergy);
    sum_average = _global_av(_index_kenergy);
    sum_ave2 = _global_av2(_index_kenergy);
    coutf << setw(12) << blk
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  // TOTAL ENERGY //////////////////////////////////////////////////////////////
  if (_measure_tenergy){
    coutf.open("../OUTPUT/total_energy.dat",ios::app);
    average  = _average(_index_tenergy);
    sum_average = _global_av(_index_tenergy);
    sum_ave2 = _global_av2(_index_tenergy);
    coutf << setw(12) << blk
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  // TEMPERATURE ///////////////////////////////////////////////////////////////
  if (_measure_temp){
    coutf.open("../OUTPUT/temperature.dat",ios::app);
    average  = _average(_index_temp);
    sum_average = _global_av(_index_temp);
    sum_ave2 = _global_av2(_index_temp);
    coutf << setw(12) << blk
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  // PRESSURE //////////////////////////////////////////////////////////////////
  if (_measure_pressure){
    coutf.open("../OUTPUT/pressure.dat",ios::app);
    average  = _average(_index_pressure);
    sum_average = _global_av(_index_pressure);
    sum_ave2 = _global_av2(_index_pressure);
    coutf << setw(12) << blk
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  // GOFR //////////////////////////////////////////////////////////////////////
  // TO BE FIXED IN EXERCISE 7
  // POFV //////////////////////////////////////////////////////////////////////
  // TO BE FIXED IN EXERCISE 4
  if (_measure_pofv){
    coutf.open("../OUTPUT/pofv.dat",ios::app);
    for (int bin = 0; bin < _n_bins_v; bin++){
      average     = _average(_index_pofv + bin) / _bin_size_v;       // media del blocco per il bin, normalizzata
      sum_average = _global_av(_index_pofv + bin) / _bin_size_v;     // somma cumulativa (normalizzata)
      sum_ave2    = _global_av2(_index_pofv + bin) / (_bin_size_v*_bin_size_v); 
  
      double v_bin = (bin + 0.5) * _bin_size_v;  // centro del bin
  
      coutf << setw(12) << blk  //indice del blocco 
      << setw(12) << v_bin        // ascissa (velocità)
      << setw(12) << average      // valore del blocco
      << setw(12) << sum_average / double(blk)  // media cumulativa
      << setw(12) << this->error(sum_average, sum_ave2, blk)    // errore statistico
      << endl; //gli ultimi due sono i valori che salvo
    }  
    coutf.close();
  }
  // MAGNETIZATION /////////////////////////////////////////////////////////////
  // TO BE FIXED IN EXERCISE 6
  // SPECIFIC HEAT /////////////////////////////////////////////////////////////
  // TO BE FIXED IN EXERCISE 6
  // SUSCEPTIBILITY ////////////////////////////////////////////////////////////
  // TO BE FIXED IN EXERCISE 6
  // ACCEPTANCE ////////////////////////////////////////////////////////////////
  double fraction;
  coutf.open("../OUTPUT/acceptance.dat",ios::app);
  if(_nattempts > 0) fraction = double(_naccepted)/double(_nattempts);
  else fraction = 0.0; 
  coutf << setw(12) << blk << setw(12) << fraction << endl;
  coutf.close();
  
  return;
}

double System :: error(double acc, double acc2, int blk){
  if(blk <= 1) return 0.0;
  else return sqrt( fabs(acc2/double(blk) - pow( acc/double(blk) ,2) )/double(blk) );
}

int System :: get_nbl(){
  return _nblocks;
}

int System :: get_nsteps(){
  return _nsteps;
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
