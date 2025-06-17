/*
  main.cpp: Coordinatore MPI per esecuzione distribuita del Genetic Algorithm sul TSP

  Questo programma gestisce l’inizializzazione e la finalizzazione dell’ambiente MPI,
  la distribuzione dei dati delle 110 città fra tutti i processi e l’esecuzione di due
  strategie di ricerca: una priva di migrazione (“_nomig”) e una con migrazione
  periodica degli individui migliori (“_mig”).

  Flusso logico:
    1. MPI_Init / MPI_Finalize: avvia e chiude l’ambiente di comunicazione parallela.
    2. Rank 0 legge da file “cap_prov_ita.dat” le coordinate delle città e ne trasmette
       il conteggio e la struttura dati a tutti gli altri processi via MPI_Bcast.
    3. Per ciascuna delle due run (disable/enable migration):
       a) setupRNG(rnd, rank): sincronizza il seme del generatore casuale per riproducibilità.
       b) Calcolo di migrate_every basato su TOTAL_GENERATIONS, TARGET_TOTAL_MIGRATIONS e size:
          definisce ogni quante generazioni invocare MPI_Allgather nel GA.
       c) Creazione dell’oggetto GA(rnd, cities) e invocazione di run(...):
          - evolveOneGeneration con selezione, crossover OX, molteplici mutazioni
            (swap, invert, shift, block‐flip, local scramble, edge optimization).
          - Se abilitata, migrazioneAllgather(n) scambia i migliori individui fra i rank.
          - Registrazione locale delle statistiche di best/avg per generazione.
       d) Raccolta su rank 0 di:
          - local_best (MPI_Gather di fitness del miglior individuo)
          - local_path (MPI_Gather del percorso relativo)
          - evolution_data_ (MPI_Gather delle statistiche per generazione)
       e) Sul rank 0:
          - Identificazione del rank con best global fitness e salvataggio del suo percorso
            in `best_path_<suffix>.dat` (chiuso ad anello).
          - Aggregazione dei vettori best_length e avg_best_half_length per generazione:
            calcolo del minimo dei best e della media degli avg su tutti i rank.
          - Scrittura in `evolution_global_<suffix>.dat` delle statistiche globali.
          - Stampa su console del riepilogo file generati.
  
  Ogni fase di comunicazione (MPI_Bcast, MPI_Gather) è progettata per minimizzare
  overhead e garantire coerenza dei dati tra i processi, mentre la struttura a due run
  permette un confronto diretto fra GA indipendenti e GA con scambio genetico.
*/

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <limits>
#include <string>
#include "random.h"
#include "ga.h"
#include "utils.h"
using namespace std;
int main(int argc, char** argv) {
    // Inizializza l'ambiente MPI
    MPI_Init(&argc, &argv);
    // Ottiene il rank (ID) del processo e il numero totale di processi
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
// Parametri di esecuzione
    const int TOTAL_GENERATIONS = 5000; // Numero totale di generazioni GA
    const int TARGET_TOTAL_MIGRATIONS = 500; // Numero desiderato di cicli di migrazione
    const int NUM_RUNS = 2;  // Esegui due run: senza e con migrazione
    bool enable_migration[NUM_RUNS] = {false, true}; // Flag per abilitare/disabilitare migrazione
    string suffixes[NUM_RUNS] = {"_nomig", "_mig"}; // Suffissi per i file di output

  // Caricamento delle città (solo il rank 0 legge da file)
    Random rnd;
    vector<City> cities;
    if (rank == 0) {
        // Il processo master carica le coordinate da file
        cities = loadCitiesFromFile("cap_prov_ita.dat");
    }
    // Comunica a tutti il numero di città lette
    int num_cities = cities.size();
    MPI_Bcast(&num_cities, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // Ridimensiona il vettore città sugli altri rank
    if (rank != 0) cities.resize(num_cities);
    MPI_Bcast(cities.data(), num_cities * sizeof(City), MPI_BYTE, 0, MPI_COMM_WORLD);
// Ciclo sui due tipi di run (no-migrazione / migrazione)
    for (int run = 0; run < NUM_RUNS; ++run) {
        // Reset RNG per mantenere condizioni iniziali identiche
        setupRNG(rnd, rank);
        
        // Calcolo migrate_every solo per run con migrazione
        int migrate_every = 0;
        if (enable_migration[run]) {
             // Distribuisce TARGET_TOTAL_MIGRATIONS equamente sui processi
            migrate_every = TOTAL_GENERATIONS / (TARGET_TOTAL_MIGRATIONS / size);
            migrate_every = max(migrate_every, 1); // Almeno ogni generazione
        }

        GA ga(rnd, cities);
        ga.run(TOTAL_GENERATIONS, migrate_every, rank, size, enable_migration[run]);

        // Raccolta best fitness locale
        double local_best = ga.getBestFitness();
        vector<double> all_bests(size);
        MPI_Gather(
            &local_best, 1, MPI_DOUBLE,  // invia il fitness locale
            all_bests.data(), 1, MPI_DOUBLE, // riceve tutti i fitness su rank 0
            0, MPI_COMM_WORLD
        );

       // --- Raccolta del percorso migliore locale da tutti i rank ---
        const auto& local_path = ga.getBestPath();
        int path_len = static_cast<int>(local_path.size());
        vector<int> all_paths;
        if (rank == 0) all_paths.resize(size * path_len);
        MPI_Gather(
            local_path.data(), path_len, MPI_INT, // invia percorso locale
            all_paths.data(), path_len, MPI_INT, // riceve tutti i percorsi su rank 0
            0, MPI_COMM_WORLD
        );

        // Salvataggio risultati globali (solo rank 0)
        if (rank == 0) {
            // Trova il rank con il miglior fitness
            int best_rank = max_element(all_bests.begin(), all_bests.end()) - all_bests.begin();

            // Estrai percorso migliore
            vector<int> best_path(path_len);
            copy_n(
                all_paths.begin() + best_rank * path_len,
                path_len,
                best_path.begin()
            );

            // Salva best_path
            string path_file = "best_path" + suffixes[run] + ".dat";
            ofstream fout(path_file);
            for (int idx : best_path) {
                fout << cities[idx].x << " " << cities[idx].y << "\n";
            }
            fout << cities[best_path[0]].x << " " << cities[best_path[0]].y << "\n";
            fout.close();
        }

        // Raccolta dati evolutivi
        auto local_evo = ga.getEvolutionData();
        int num_generations = static_cast<int>(local_evo.size());

        // Preparazione buffer locali
        vector<int> local_gens(num_generations);
        vector<double> local_best_lens(num_generations), local_avg_lens(num_generations);
        for (int g = 0; g < num_generations; ++g) {
            local_gens[g] = local_evo[g].generation;
            local_best_lens[g] = local_evo[g].best_length;
            local_avg_lens[g] = local_evo[g].avg_best_half_length;
        }

        // Raccolta dati globali
        vector<int> all_gens;
        vector<double> all_bests_evo, all_avgs;
        if (rank == 0) {
            all_gens.resize(num_generations * size);
            all_bests_evo.resize(num_generations * size);
            all_avgs.resize(num_generations * size);
        }
 // Gather generation indices 
        MPI_Gather(
            local_gens.data(), num_generations, MPI_INT,
            rank == 0 ? all_gens.data() : nullptr, num_generations, MPI_INT,
            0, MPI_COMM_WORLD
        );
         // Gather best_length per generazione

        MPI_Gather(
            local_best_lens.data(), num_generations, MPI_DOUBLE,
            rank == 0 ? all_bests_evo.data() : nullptr, num_generations, MPI_DOUBLE,
            0, MPI_COMM_WORLD
        );
// Gather avg_best_half_length per generazione
        MPI_Gather(
            local_avg_lens.data(), num_generations, MPI_DOUBLE,
            rank == 0 ? all_avgs.data() : nullptr, num_generations, MPI_DOUBLE,
            0, MPI_COMM_WORLD
        );

        // Salvataggio evolution_global (solo rank 0)
        if (rank == 0) {
            string evo_file = "evolution_global" + suffixes[run] + ".dat";
            ofstream fevo(evo_file);
            fevo << "# Generation Best_Length AvgBestHalf\n";

            for (int g = 0; g < num_generations; ++g) {
                double global_best = numeric_limits<double>::infinity();
                double sum_avg = 0.0;
                 // Per ogni rank, aggiorna il best globale e somma gli avg
                for (int r = 0; r < size; ++r) {
                    int idx = r * num_generations + g;
                    global_best = min(global_best, all_bests_evo[idx]);
                    sum_avg += all_avgs[idx];
                }
                // Scrive generazione, best globale e media degli avg
                fevo << g << " " << global_best << " " << (sum_avg / size) << "\n";
            }
            fevo.close();
  // Solo rank 0 stampa a video i file prodotti
            cout << "Run " << suffixes[run] << " completato!\n";
            cout << "- Percorso migliore: best_path" << suffixes[run] << ".dat\n";
            cout << "- Dati evolutivi: evolution_global" << suffixes[run] << ".dat\n";
        }
    }

    MPI_Finalize();
    return 0;
}