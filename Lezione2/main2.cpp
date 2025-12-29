#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "random.h"

using namespace std;

/* ================= FUNZIONE PER ERRORE =====================
   Restituisce l'errore statistico tramite blocking method:
   σ = sqrt( (<x²> – <x>²) / N )
============================================================= */
double error(const vector<double>& AV, const vector<double>& AV2, int n) {
    if (n == 0) return 0.0;
    return sqrt((AV2[n] - AV[n] * AV[n]) / n);
}



int main() {

    /***********************************************************
     *    INIZIALIZZAZIONE GENERATORE RANDOM
     ***********************************************************/
    Random rnd;
    int seed[4], p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()) Primes >> p1 >> p2;
    else cerr << "ERRORE apertura Primes\n";
    Primes.close();

    ifstream input("seed.in");
    string property;
    if (input.is_open()) {
        while (!input.eof()) {
            input >> property;
            if (property == "RANDOMSEED")
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
        }
        rnd.SetRandom(seed, p1, p2);
        input.close();
    }
    else cerr << "ERRORE apertura seed.in\n";



/***************************************************************
RANDOM WALK LATTICE 3D (DISCRETO) 
 Calcolo di ⟨r²(n)⟩ e √⟨r²(n)⟩ con blocking
***************************************************************/
    int M = 10000;           // # Random Walk totali
    int N = 100;             // # Blocchi
    int L = M/N;             // # RW per blocco
    int n_step = 100;        // lunghezza delle passeggiate

    vector<double> ave(n_step), ave2(n_step), err_prog(n_step);


/*** SIMULAZIONE A BLOCCHI ***/
    for(int i=0; i<N; i++){

        vector<double> r2_block(n_step, 0.0);

        for(int j=0; j<L; j++){ // L RW per blocco

            vector<int> pos = {0,0,0}; // → riparto sempre dall’origine

            for(int k=0; k<n_step; k++){

                int direzione = rnd.Rannyu(0,3);   // 0=x 1=y 2=z
                int verso = rnd.Rannyu(0,2);       // +1 o -1

                if(verso==0) pos[direzione]--;
                else         pos[direzione]++;

                r2_block[k] += pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2];
            }
        }

        for(int k=0; k<n_step; k++){
            r2_block[k] /= L;
            ave[k]  += r2_block[k];
            ave2[k] += r2_block[k]*r2_block[k];
        }
    }


/*** Calcolo statistiche finali ***/
    for(int k=0; k<n_step; k++){
        ave[k]  /= N;
        ave2[k] /= N;
        err_prog[k] = sqrt((ave2[k]-ave[k]*ave[k])/(N-1));
    }


/*** Salvataggio risultati r² discreto ***/
    ofstream out("RW_discreto.txt");
    for(int i=0; i<n_step; i++){
        double r = sqrt(ave[i]);
        double err = 0.5 * err_prog[i] / r;
        out << i+1 << "\t" << r << "\t" << err << endl;
    }
    out.close();



/*** Salvataggio POSIZIONI PER PLOT 3D ***/
    ofstream outplot1("output_disegnoRW.txt");
    vector<int> posizione = {0,0,0}; // random walk singolo SOLO PER GRAFICO

    for (int k=0; k<n_step; k++){
        int d = rnd.Rannyu(0,3);
        int v = rnd.Rannyu(0,2);
        if(v==0) posizione[d]--; else posizione[d]++;
        outplot1 << posizione[0] << "\t" << posizione[1] << "\t" << posizione[2] << endl;
    }
    outplot1.close();




/***************************************************************
 RANDOM WALK 3D CONTINUO
 Passi di modulo 1 all’angolo casuale
***************************************************************/
    vector<double> ave_c(n_step), ave2_c(n_step), err_prog_c(n_step);


/*** Simulazione continua con direzione uniforme sulla sfera ***/
    for(int i=0; i<N; i++){

        vector<double> r2_block(n_step,0.0);

        for(int j=0; j<L; j++){

            vector<double> pos = {0.,0.,0.};

            for(int k=0; k<n_step; k++){

                double theta = rnd.Rannyu(0,M_PI);
                double phi   = rnd.Rannyu(0,2*M_PI);

                pos[0] += sin(theta)*cos(phi);
                pos[1] += sin(theta)*sin(phi);
                pos[2] += cos(theta);

                r2_block[k] += pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2];
            }
        }

        for(int k=0; k<n_step; k++){
            r2_block[k] /= L;
            ave_c[k]  += r2_block[k];
            ave2_c[k] += r2_block[k]*r2_block[k];
        }
    }


/*** Analisi statistica continuo ***/
    for(int k=0; k<n_step; k++){
        ave_c[k]  /= N;
        ave2_c[k] /= N;
        err_prog_c[k] = sqrt((ave2_c[k]-ave_c[k]*ave_c[k])/(N-1));
    }


/*** Output dati finali ***/
    ofstream outc("RW_continuo.txt");
    for(int i=0; i<n_step; i++){
        double r = sqrt(ave_c[i]);
        double err = 0.5 * err_prog_c[i] / r;
        outc << i+1 << "\t" << r << "\t" << err << endl;
    }
    outc.close();



/*** RW per grafico 3D continuo ***/
    ofstream outplot2("output_disegnoRW_c.txt");
    vector<double> posizione_c = {0.,0.,0.};

    for(int k=0; k<n_step; k++){
        double theta = rnd.Rannyu(0,M_PI);
        double phi   = rnd.Rannyu(0,2*M_PI);
        posizione_c[0]+=sin(theta)*cos(phi);
        posizione_c[1]+=sin(theta)*sin(phi);
        posizione_c[2]+=cos(theta);
        outplot2 << posizione_c[0] << "\t" << posizione_c[1] << "\t" << posizione_c[2] << endl;
    }
    outplot2.close();


    rnd.SaveSeed();
    return 0;
}
