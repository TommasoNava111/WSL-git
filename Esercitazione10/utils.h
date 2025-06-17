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
#include <vector>


struct City;
std::vector<City> loadCitiesFromFile(const std::string& filename);
// Legge due numeri primi da file per inizializzare il generatore.
// Scorre fino alla riga corrispondente al `rank`.
void ReadPrimes(const std::string& filename, int& p1, int& p2, int rank);

// Estrae i semi da seedFile e chiama rnd.SetRandom(...)
void InitializeRandom(Random& rnd, const std::string& seedFile, int p1, int p2);

// Configura il generatore globale per un dato `rank`
// (chiama ReadPrimes(filename, p1,p2, rank) e InitializeRandom).
void setupRNG(Random& rnd, int rank);

#endif // UTILS_H