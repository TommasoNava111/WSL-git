// ============================================================================
// utils.h
//
// Header che definisce funzioni di utilità per la gestione del generatore
// di numeri casuali. Serve per inizializzare correttamente il generatore
// (oggetto della classe Random), leggendo semi e numeri primi da file.
//
// Utilizzato nel contesto dell'esercizio 09.1 per inizializzare il generatore
// `rnd` usato durante il Genetic Algorithm (GA) per il TSP.
//
// Le funzioni contenute qui permettono di rendere l'inizializzazione portabile,
// modulare e leggibile.
//
// ============================================================================
#ifndef UTILS_H
#define UTILS_H

#include <string>
#include "random.h"

// Legge due numeri primi da file per inizializzare il generatore
void ReadPrimes(const std::string& filename, int& p1, int& p2);

// Estrae i semi da seedFile e chiama rnd.SetRandom(...)
void InitializeRandom(Random& rnd, const std::string& seedFile, int p1, int p2);

// Configura il generatore globale (chiama ReadPrimes e InitializeRandom)
void setupRNG(Random& rnd);

#endif // UTILS_H
