// ============================================================================
// ga.h
//
// Header per la classe GA che implementa un algoritmo genetico per il TSP.
// Contiene:
// - Dichiarazione della struttura City (coordinate x,y).
// - Parametri fissi del GA (numero di città, dimensione popolazione, ecc.).
// - Prototipi delle funzioni per generazione, selezione, mutazioni,
//   crossover, 2-opt e valutazione fitness.
// ============================================================================

#ifndef GA_H
#define GA_H

#include <vector>
#include <string>
#include "random.h"

using namespace std;  
// Struttura per salvare le coordinate di una città
struct City { 
    double x, y; 
};

class GA {
public:
    // Costruttore: riceve il generatore casuale e la flag on_circle
    GA(Random& rnd, bool on_circle);
    // Esegue l’intero ciclo del GA e salva i risultati su file
    void run();
    
    // Operatore di mutazione: shift circolare
    void shiftMut(vector<int>& p);
private:
    // Genera le coordinate delle N_CITIES (su cerchio o quadrato)
    void generateCities();
    // Inizializza la popolazione con permutazioni casuali valide
    void initializePopulation();
    // Esegue una sola generazione (selezione, crossover, mutazioni, elitismo)
    void evolveOneGeneration();
    // Calcola la fitness (1/lunghezza del percorso)
    double calculateFitness(const vector<int>& path) const;
    // Verifica che un percorso sia valido (ogni città una sola volta)
    bool isValid(const vector<int>& p) const;
    // Crossover Order-Xpreserving tra due genitori
    vector<int> crossoverOX(const vector<int>& A, const vector<int>& B);
    // Mutazioni base: swap e inversione
    void swapMut(vector<int>& p);
    void invertMut(vector<int>& p);
    // Distanza tra due città (L2 o L1 a seconda di USE_L2)
    double dist(int a, int b) const;
    // Mutazione block swap
    void blockSwapMut(vector<int>& p);
    double selection_p_;  // Esponente di selezione
    vector<size_t> sorted_indices_;  // Indici ordinati per fitness
     // Ordina popolazione per fitness
    void orderPopulation();
    // Mutazione local scramble 
    void localScrambleMut(vector<int>& p);

    // Parametri costanti del GA
    static constexpr int    N_CITIES      = 34;  // Numero di città
    static constexpr int    POP_SIZE      = 500; // Dimensione popolazione
    static constexpr double MUT_RATE      = 0.08;  // Probabilità mutazione
    static constexpr double CROSS_RATE    = 0.7; // Probabilità crossover
    static constexpr int    GENERATIONS   = 5000;   // Numero di generazioni
    static constexpr bool   USE_L2        = true;   // Usa L2 (altrimenti L1)

    // Stato interno
    Random&                    rnd_;          // Generatore di numeri casuali
    bool                       on_circle_;    // true → città su circonferenza
    string                     suffix_;       // suffisso per file di output
    vector<City>               cities_;       // coordinate delle città
    vector<vector<int>>        population_;   // percorsi
    vector<double>             fitness_;      // fitness di ciascun individuo
};

#endif // GA_H
