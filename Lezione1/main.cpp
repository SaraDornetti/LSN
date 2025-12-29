#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "random.h"
using namespace std;

/************* FUNZIONI *************/

// Calcola errore statistico (incertezza sulla media)
double error(const vector<double>& AV, const vector<double>& AV2, int n){
    if(n==0) return 0.;  // Se non ci sono dati, errore = 0
    return sqrt((AV2[n]-AV[n]*AV[n])/n);  
}

// Salva 3 colonne su file: blocco || media || errore
void SaveData(const string &filename, const vector<double>& ave, const vector<double>& err){
    ofstream out(filename);          // apro file in scrittura
    if(!out){ cerr<<"Errore apertura file "<<filename<<endl; exit(1);}
    for(size_t i=0;i<ave.size();i++) 
        out<<i+1<<"\t"<<ave[i]<<"\t"<<err[i]<<endl;  // scrivo indice blocco, media, errore
}

// Esegue medie a blocchi
void BlockAverage(Random& rnd,int M,int N,const string& file_mean,const string& file_variance){

   // M = numero di lanci
   // N = numero di blocchi
   int L = M/N;  // numero di lanci per blocco

    //inizializzo a zero vettori di N componenti per contenere medie, medie dei quadrati, somme ed errore
    vector<double> ave(N,0.), ave2(N,0.), sum(N,0.), sum2(N,0.), err(N,0.); // per media uniforme
    vector<double> aveV(N,0.), ave2V(N,0.), sumV(N,0.), sum2V(N,0.), errV(N,0.); // per varianza

    for(int i=0;i<N;i++){                 // ciclo sui blocchi
        double s1=0, s2=0;                // variabili per calcolare le somme inizializzate a zero
        for(int j=0;j<L;j++){             // ciclo sul numero di lanci per blocco
            double r = rnd.Rannyu();      // estrazione uniforme
            s1 += r;                      // accumulo media
            s2 += pow(r-0.5,2);           // accumulo varianza (centrata)
        }
        ave[i]  = s1/L;             // media
        ave2[i]  = ave[i]*ave[i];   // media al quadrato

        aveV[i] = s2/L;             // media varianza 
        ave2V[i] = aveV[i]*aveV[i]; // media varianza al quadrato
    }

    // calcolo medie progressive + errori
    for(int i=0;i<N;i++){
        for(int j=0;j<=i;j++){
            sum[i]+=ave[j];         // somme medie progressive
            sum2[i]+=ave2[j];       // somme medie al quadrato progressive

            sumV[i]+=aveV[j];       // somme medie varianza progressive
            sum2V[i]+=ave2V[j];     // somme medie al quadrato varianza progressive
        }
      sum[i]/=(i+1);
      sum2[i]/=(i+1); 
      err[i]=error(sum,sum2,i);
      
      sumV[i]/=(i+1); 
      sum2V[i]/=(i+1); 
      errV[i]=error(sumV,sum2V,i);
    }

    SaveData(file_mean,sum,err);           // salvo media uniforme
    SaveData(file_variance,sumV,errV);     // salvo media varianza
}

// Test χ² su distribuzione uniforme
void ChiSquareTest(Random& rnd,int lanci,int bins,const string& file){
    ofstream out(file);
    double atteso=(double)lanci/bins;    // valore atteso in ogni bin
    int N_att=500;                           // numero di 𝝌² tests da fare   
    for(int t=0;t<N_att;t++){                  // ripeto il test 100 volte
        vector<int> count(bins,0);           // contatore per ogni bin
        for(int i=0;i<lanci;i++) 
            count[int(rnd.Rannyu()*bins)]++; // estrazione e incremento bin
        double chi2=0;
        for(int i=0;i<bins;i++) 
            chi2+=pow(count[i]-atteso,2)/atteso; // formula χ²
        out<<t+1<<"\t"<<chi2<<endl;           // salvo risultato
    }
}

// Central Limit Theorem per diversi valori di N
void CLT(Random& rnd,int M,const vector<int>& valoriN,const string& base){
    for(int N: valoriN){
        ofstream out(base+"_N"+to_string(N)+".txt"); // nome file dinamico
        for(int i=0;i<M;i++){
            double sU=0,sE=0,sC=0;
            for(int j=0;j<N;j++){
                sU+=rnd.Rannyu();          // somma numeri uniformi
                sE+=rnd.Exponential(1.);   // somma numeri esponenziali
                sC+=rnd.CauchyLorentz(1.,0.); // somma numeri Cauchy-Lorentz
            }
            out<<sU/N<<"\t"<<sE/N<<"\t"<<sC/N<<"\n";  // media su N estrazioni
        }
    }
}

// Buffon Needle esteso su più linee
void Buffon(Random& rnd,int M,int N,double L_ago,double d,int nLines,const string& file){
    int L=M/N;                               // lanci per blocco
    vector<double> ave(N,0.), ave2(N,0.), sum(N,0.), sum2(N,0.), err(N,0.);

    // ciclo sui blocchi
    for(int i=0;i<N;i++){
        int hit=0;                            // numero intersezioni
        for(int j=0;j<L;j++){
            double y=rnd.Rannyu(0,d*nLines);  // centro dell'ago
            int k=y/d;                        // riga in cui cade
            double a=rnd.Angle();             // angolo
            if( (y+(L_ago/2)*sin(a)>d*(k+1)) || (y-(L_ago/2)*sin(a)<d*k)) hit++; // intersezione
        }
        ave[i]=2*L_ago*L/(hit*d);    
        ave2[i]=ave[i]*ave[i];
    }

    // medie progressive + errore
    for(int i=0;i<N;i++){
        for(int j=0;j<=i;j++){ sum[i]+=ave[j]; sum2[i]+=ave2[j]; }
        sum[i]/=(i+1); sum2[i]/=(i+1); err[i]=error(sum,sum2,i);
    }
    SaveData(file,sum,err);                  // salvo risultati
}


/******************* MAIN *******************/
int main(int argc,char* argv[]){
    
    Random rnd;           
    int seed[4],p1,p2;

    // apro file Primes per inizializzazione RNG
    ifstream Primes("Primes"); Primes>>p1>>p2; Primes.close();

    // apro file seed.in per inizializzare RNG
    ifstream input("seed.in"); string prop;
    while(input>>prop) 
        if(prop=="RANDOMSEED") input>>seed[0]>>seed[1]>>seed[2]>>seed[3];
    input.close(); 
    rnd.SetRandom(seed,p1,p2);  // inizializzo RNG

    // ===========================================
    // esercizio 1.1.1 e 1.1.2: media e varianza
    // ===========================================
    BlockAverage(rnd,10000,100,"output.txt","output_v.txt");

    // ===========================================
    // esercizio 1.1.3: test chi quadro
    // ===========================================
    ChiSquareTest(rnd,10000,100,"output_ChiQuadro.txt");

    // ===========================================
    // esercizio 1.2: Teorema del limite centrale
    // ===========================================
    CLT(rnd,10000,{1,2,10,100},"output_S");

    // ===========================================
    // esercizio 1.3: Esperimento Buffon
    // ===========================================
    Buffon(rnd,100000,100,0.34,1.,10,"output_Buffon.txt");

    rnd.SaveSeed();  
    return 0;
}
