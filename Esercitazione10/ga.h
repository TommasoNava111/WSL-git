#ifndef GA_H
#define GA_H
/*
  Genetic Algorithm per il TSP con supporto MPI

  Contiene la definizione della classe GA e delle strutture dati associate:
    - City: rappresenta le coordinate (x,y) di una città.
    - GenerationData: raccoglie statistiche per generazione (numero, lunghezza migliore, lunghezza media della metà migliore).
  La classe GA offre:
    - Costruttore che riceve un generatore di numeri casuali e la lista di città.
    - run(): esegue l’evoluzione per un dato numero di generazioni, con facoltativa migrazione dei migliori individui fra più processi MPI.
    - Metodi per ottenere il percorso e il fitness migliore trovati.
    - Member functions per inizializzare la popolazione, evolvere (selezione, crossover, mutazioni), migrare individui e calcolare fitness.
    - Vari operatori di mutazione e crossover OX per il TSP.
    - Pre-calcolo e accesso alla matrice delle distanze Euclidee fra città.
  
  Parametri di configurazione (numero di città, dimensione popolazione, tassi di mutazione/crossover, selezione) sono definiti come constexpr.
*/

#include <vector>
#include <string>
#include "random.h"
#include <mpi.h>
using namespace std;
// Struttura per contenere le coordinate di una città
struct City { 
    double x, y; 
};
// Struttura per memorizzare i dati di una generazione (per output)
struct GenerationData {
    int generation;
    double best_length;
    double avg_best_half_length;
};
class GA {
public:
    GA(Random& rnd, const vector<City>& cities);

    // Esegue l’evoluzione per total_generations con migrazione ogni migrate_every
    // rank e size servono per MPI
     void run(int total_generations, int migrate_every, int rank, int size, bool enable_migration);

    // Recupera il migliore fitness e percorso trovati da questo rank
    double getBestFitness() const;
    const vector<int>& getBestPath() const;
    // Vettore pubblico per recuperare i dati di evoluzione
    vector<GenerationData> evolution_data_;
    vector<GenerationData> getEvolutionData() const { return evolution_data_; }

private:
    // Inizializza popolazione e fitness
    void initializePopulation();

    // Un singolo passo di evoluzione (selezione, crossover, mutazioni)
    void evolveOneGeneration();

    // Migrazione collettiva tramite MPI_Allgather
    void migrationAllgather(int n_people);

    // Utility per ordinare e calcolare fitness
    void orderPopulation();
    double calculateFitness(const vector<int>& path) const;

    // Controllo validità e operatori GA
    bool isValid(const vector<int>& p) const;
    vector<int> crossoverOX(const vector<int>& A,
                                 const vector<int>& B);
  // Operatori di mutazione                      
    void swapMut(vector<int>& p);
    void invertMut(vector<int>& p);
    void shiftMut(vector<int>& p);
    void blockSwapMut(vector<int>& p);
    void localScrambleMut(vector<int>& p);
    void edgeOptimizerMut(vector<int>& p);
    void blockFlipMut(vector<int>& p);
    vector<vector<double>> distance_matrix_;
    
    void precomputeDistances();   // calcola matrice delle distanze
    vector<int> getRandomIndices(int n);
   

    // Calcolo distanza al quadrato (L2)
    double dist(int a, int b) const;

    // Restituisce gli indici dei migliori (take_best=true) o peggiori (false)
   vector<int> getTopIndices(int n, bool take_best = true);

    // Parametri del GA
    static constexpr int    N_CITIES    = 110;
    static constexpr int    POP_SIZE    = 500;
    static constexpr double MUT_RATE    = 0.08;
    static constexpr double CROSS_RATE  = 0.8;
    static constexpr double SELECTION_P = 5.0;

    // Stato interno
    Random& rnd_;
    vector<City> cities_;
    vector<vector<int>> population_;
    vector<double> fitness_;
    vector<size_t> sorted_indices_;
    int rank_, size_;

    // ** Record locale del migliore **
    double best_fitness_;
    vector<int> best_path_;
};

#endif // GA_H