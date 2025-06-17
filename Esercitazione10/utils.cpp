// ============================================================================
// utils.cpp
//
// Implementazione delle funzioni dichiarate in utils.h per l'inizializzazione
// del generatore di numeri casuali. Utilizza due file esterni:
// - "Primes": contenente due numeri primi per l'inizializzazione
// - "seed.in": contenente i semi iniziali per il generatore
//
// Queste funzioni sono fondamentali per ottenere un comportamento riproducibile
// e casuale nella simulazione GA.
//
// ============================================================================

#include "utils.h"
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <vector>
#include "ga.h"

using namespace std;

// Funzione che legge due numeri primi da un file, necessari per rnd.SetRandom
void ReadPrimes(const string& filename, int& p1, int& p2, int rank) {
    ifstream filePrimes(filename);
    if (!filePrimes.is_open()) {
        cerr << "ERRORE: Impossibile aprire " << filename << "\n";
        exit(1);
    }
    // Scorre fino alla riga corrispondente al rank
    for (int i = 0; i <= rank; ++i) {
        filePrimes >> p1 >> p2;
    }
    filePrimes.close();
}
void setupRNG(Random& rnd, int rank) {
    int p1, p2;
    ReadPrimes("Primes", p1, p2, rank);
    InitializeRandom(rnd, "seed.in", p1, p2);
}

// Inizializza l'oggetto rnd leggendo 4 semi interi da un file e chiamando SetRandom
void InitializeRandom(Random& rnd, const string& seedFile, int p1, int p2) {
    ifstream fileSeed(seedFile);
    if (!fileSeed.is_open()) {
        cerr << "ERRORE: Impossibile aprire " << seedFile << "\n";
        exit(1);
    }

    string keyword;
    int seedArray[4]; // array che conterrà i 4 semi

    // Cerca la parola chiave "RANDOMSEED" e legge i valori seguenti
    while (fileSeed >> keyword) {
        if (keyword == "RANDOMSEED") {
            fileSeed 
              >> seedArray[0]
              >> seedArray[1]
              >> seedArray[2]
              >> seedArray[3];
            rnd.SetRandom(seedArray, p1, p2); // inizializza il generatore
            break;
        }
    }
    fileSeed.close();
}