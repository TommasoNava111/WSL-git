// Es02.2.cpp
//
// Exercise 02.2: 3D Random Walks on a cubic lattice and in the continuum
// --------------------------------------------------------------------------------
// Scopo: simulare cammini casuali (RW) 3D partendo dall'origine, ripetuti M volte
// (organizzati in N_BLOCKS blocchi di K_WALKS cammini ciascuno), per verificare il
// comportamento diffusivo e calcolare la dipendenza di √⟨|r_N|²⟩ dal numero di passi N.
// 
// 1) Reticolo cubico: passo di lunghezza a=1 lungo x,y o z con direzione +/- scelta
//    uniformemente. Ad ogni passo si aggiorna la posizione (x,y,z) e si calcola r².
// 
// 2) Cammino continuo: passo di lunghezza a=1 in direzione casuale ottenuta tramite
//    θ ∈ [0,π], φ ∈ [0,2π] campionati uniformemente sul solido (θ da inversione 1−2u,
//    φ=2πv). Si accumula r² ad ogni passo.
// 
// Metodo:
// - Dividere M=TOT_WALKS=N_BLOCKS*K_WALKS estrazioni in blocchi di K_WALKS cammini.
// - Per ciascun blocco e per ogni passo i=0..N_STEPS-1 calcolare la media di r² su
//   K_WALKS cammini (block mean), e accumulare somme e somme quadratiche per i blocchi.
// - Al termine, per ogni passo i, ottenere la media cumulativa ⟨r²⟩ e l’errore della media
//   (deviazione standard della media), trasformando poi in √⟨r²⟩ con propagazione errore.
// - Salvare i risultati in lattice.txt e continuum.txt: righe [step] [√⟨r²⟩] [errore].
// 
// Risultato atteso: √⟨|r_N|²⟩ ∝ √N (diffusione libera), verificabile tramite fit k√N.

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include "random.h"
using namespace std;

typedef double Real;

// Parametri di simulazione
const int N_BLOCKS = 100;   // numero di blocchi
const int K_WALKS  = 100;   // cammini per blocco
const int N_STEPS  = 100;   // passi per cammino

// Prototipi funzioni
void ReadPrimes(const string& filename, int& p1, int& p2);
void InitializeRandom(Random& rnd, const string& seedFile, int p1, int p2);
void ComputeRandomWalkStats(Random& rnd,
                            int nBlocks,
                            int kWalks,
                            const string& outFileLat,
                            const string& outFileCont);

int main() {
    Random rnd;
    int p1, p2;
    ReadPrimes("Primes", p1, p2);
    InitializeRandom(rnd, "seed.in", p1, p2);

    ComputeRandomWalkStats(rnd, N_BLOCKS, K_WALKS, "Es02.2ret.txt", "Es02.2con.txt");

    rnd.SaveSeed();
    return 0;
}

// Funzione che esegue il calcolo a blocchi per i random walk
void ComputeRandomWalkStats(Random& rnd,
                            int nBlocks,
                            int kWalks,
                            const string& outFileLat,
                            const string& outFileCont) {
    // Array cumulativi per r^2
    Real cumSumLat[N_STEPS] = {0.0}, cumSum2Lat[N_STEPS] = {0.0};
    Real cumSumCont[N_STEPS] = {0.0}, cumSum2Cont[N_STEPS] = {0.0};

    // Ciclo sui blocchi
    for (int block = 0; block < nBlocks; ++block) {
        // Sommatorie di blocco per r^2
        Real blockSumLat[N_STEPS] = {0.0};
        Real blockSumCont[N_STEPS] = {0.0};

        // Simulazione di kWalks cammini
        for (int walk = 0; walk < kWalks; ++walk) {
            // Inizio cammino
            Real x_lat=0, y_lat=0, z_lat=0;
            Real x_cont=0, y_cont=0, z_cont=0;
            blockSumLat[0] += 0.0;
            blockSumCont[0] += 0.0;

            for (int step = 1; step < N_STEPS; ++step) {
                // Reticolo cubico
                int axis = int(rnd.Rannyu() * 3);
                int sign = (rnd.Rannyu() < 0.5) ? 1 : -1;
                if      (axis == 0) x_lat += sign;
                else if (axis == 1) y_lat += sign;
                else                z_lat += sign;
                Real r2_lat = x_lat*x_lat + y_lat*y_lat + z_lat*z_lat;
                blockSumLat[step] += r2_lat;

                // Passo continuo
                Real u = rnd.Rannyu();
                Real v = rnd.Rannyu();
                Real theta = acos(1.0 - 2.0*u);
                Real phi   = 2.0 * M_PI * v;
                Real dx = sin(theta) * cos(phi);
                Real dy = sin(theta) * sin(phi);
                Real dz = cos(theta);
                x_cont += dx; y_cont += dy; z_cont += dz;
                Real r2_cont = x_cont*x_cont + y_cont*y_cont + z_cont*z_cont;
                blockSumCont[step] += r2_cont;
            }
        }

        // Aggiornamento cumulativi con medie di blocco
        for (int step = 0; step < N_STEPS; ++step) {
            Real meanLat = blockSumLat[step] / kWalks;
            cumSumLat[step]  += meanLat;
            cumSum2Lat[step] += meanLat * meanLat;
            Real meanCont = blockSumCont[step] / kWalks;
            cumSumCont[step]  += meanCont;
            cumSum2Cont[step] += meanCont * meanCont;
        }
    }

    // Apertura file di output
    ofstream outLat(outFileLat.c_str());
    ofstream outCont(outFileCont.c_str());
    if (!outLat || !outCont) {
        cerr << "ERRORE: impossibile aprire file di output" << endl;
        exit(1);
    }

    // Calcolo medie finali e incertezze, scrittura su file
    for (int step = 0; step < N_STEPS; ++step) {
        // Reticolo cubico
        Real m2Lat = cumSumLat[step] / nBlocks;
        Real var2Lat = (cumSum2Lat[step]/nBlocks - m2Lat*m2Lat) / (nBlocks - 1);
        Real err2Lat = sqrt(var2Lat) / sqrt(nBlocks);
        Real mLat = sqrt(m2Lat);
        Real errLat = (mLat > 0) ? err2Lat / (2.0 * mLat) : 0.0;

        // Continuo
        Real m2Cont = cumSumCont[step] / nBlocks;
        Real var2Cont = (cumSum2Cont[step]/nBlocks - m2Cont*m2Cont) / (nBlocks - 1);
        Real err2Cont = sqrt(var2Cont) / sqrt(nBlocks);
        Real mCont = sqrt(m2Cont);
        Real errCont = (mCont > 0) ? err2Cont / (2.0 * mCont) : 0.0;

        outLat  << step << " " << mLat  << " " << errLat  << "\n";
        outCont << step << " " << mCont << " " << errCont << "\n";
    }
}

// Lettura primi per il RNG
void ReadPrimes(const string& filename, int& p1, int& p2) {
    ifstream file(filename.c_str());
    if (!file) { cerr << "ERRORE: Primes non trovato" << endl; exit(1); }
    file >> p1 >> p2;
    file.close();
}

// Inizializzazione RNG con seed
void InitializeRandom(Random& rnd, const string& seedFile, int p1, int p2) {
    ifstream file(seedFile.c_str());
    if (!file) { cerr << "ERRORE: seed.in non trovato" << endl; exit(1); }
    string key;
    int seedArr[4];
    while (file >> key) {
        if (key == "RANDOMSEED") {
            file >> seedArr[0] >> seedArr[1] >> seedArr[2] >> seedArr[3];
            rnd.SetRandom(seedArr, p1, p2);
            break;
        }
    }
    file.close();
}
