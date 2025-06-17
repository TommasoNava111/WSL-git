// ============================================================================
// main.cpp
//
// Punto di ingresso del programma per il TSP con Genetic Algorithm.
// Flusso principale:
//  1. Inizializza RNG.
//  2. Esegue GA per 34 città su circonferenza.
//  3. Esegue GA per 34 città in quadrato.
//  4. Salva stato finale del generatore per riproducibilità.
// ============================================================================

#include <iostream>
#include "random.h"
#include "utils.h"
#include "ga.h"

using namespace std;  // Per cout, endl

int main() {
    // 1) Configurazione del generatore di numeri casuali
    Random rnd;
    setupRNG(rnd);  // legge Primes e seed.in, chiama rnd.SetRandom

    // 2) Esecuzione GA per disposizione circolare delle città
    {
        GA ga_circle(rnd, true /*on_circle*/);
        ga_circle.run();
        cout << "Simulazione su cerchio completata." << endl;
    }

    // 3) Esecuzione GA per disposizione quadrata delle città
    {
        GA ga_square(rnd, false /*on_circle*/);
        ga_square.run();
        cout << "Simulazione su quadrato completata." << endl;
    }

    // 4) Salvataggio finale dei semi del RNG per riproducibilità
    rnd.SaveSeed();

    return 0;
}
