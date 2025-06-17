/*
  Implementazione del Genetic Algorithm e del dialogo con MPI

  Questo file definisce tutti i metodi dichiarati in ga.h:
    - precomputeDistances(): costruisce la matrice NxN delle distanze Euclidee.
    - initializePopulation(): genera POP_SIZE percorsi casuali validi e calcola il fitness iniziale.
    - run(): ciclo principale con evolveOneGeneration(), calcolo delle statistiche di best e avg, migrazione collettiva (MPI_Allgather) e scrittura su file dei risultati (evoluzione e percorso migliore).
    - evolveOneGeneration(): selezione esponenziale, crossover OX, più mutazioni (swap, invert, shift, block‐flip, local scramble, ottimizzazione di bordo), con elitismo.
    - migrationAllgather(): scambio periodico dei migliori n individui fra processi MPI, sostituendo i peggiori con nuovi immigrati.
    - Funzioni di utilità: orderPopulation(), calculateFitness(), getTopIndices(), isValid(), crossoverOX(), mutazioni varie (swapMut, invertMut, shiftMut, blockSwapMut, localScrambleMut, edgeOptimizerMut, blockFlipMut).
    - loadCitiesFromFile(): helper per caricare da file le coordinate delle città.
  
  La logica principale è focalizzata su un GA parallelo “continentale” in cui ogni rank MPI evolve una popolazione indipendente e, periodicamente, migra i migliori individui per favorire la diversità genetica e accelerare la convergenza.
*/

#include "ga.h"
#include <algorithm>
#include <numeric>
#include <fstream>
#include <cmath>
#include <mpi.h>
#include"utils.h"
#include <vector>
#include <string>
using namespace std;
// Costruttore: riceve un generatore random e la lista delle città
GA::GA(Random& rnd, const vector<City>& cities)
    : rnd_(rnd), cities_(cities)  // inizializza i membri rnd_ e cities_
{
    precomputeDistances();  // riempie la matrice delle distanze fra tutte le coppie di città
    initializePopulation();   // genera POP_SIZE percorsi casuali e ne calcola il fitness
    // inizializza il best locale
    int ib = max_element(fitness_.begin(), fitness_.end()) - fitness_.begin();
    best_fitness_ = fitness_[ib]; // conserva il valore di fitness
    best_path_    = population_[ib]; // conserva il percorso corrispondente
} 

// Metodo principale: esegue l’evoluzione su total_generations generazioni
// migrate_every = intervallo di generazioni dopo cui migrare, rank e size per MPI
// enable_migration = abilita/disabilita la fase di migrazione
void GA::run(int total_generations, int migrate_every, int rank, int size, bool enable_migration) {
    rank_ = rank;  // memorizza il proprio rank MPI
    size_ = size;  // memorizza il numero di processi MPI
   string suffix = enable_migration ? "_mig" : "_nomig"; // suffisso per nomi file

    // Prepara il vettore delle statistiche
    evolution_data_.clear();
    evolution_data_.reserve(total_generations);

    // 1) Ciclo evolutivo
    for (int gen = 0; gen < total_generations; ++gen) {
        evolveOneGeneration(); // un passo di selezione, crossover e mutazioni

        // Ordina e calcola best e avg della metà migliore
        orderPopulation();
        double best_len = 1.0 / fitness_[sorted_indices_[0]];
        double sum = 0.0;
        for (int i = 0; i < POP_SIZE/2; ++i) {
            sum += 1.0 / fitness_[sorted_indices_[i]];
        }
        double avg_best_half = sum / (POP_SIZE/2);

        // Memorizza
        evolution_data_.push_back({gen, best_len, avg_best_half});

        // Se abilitata migrazione e arrivati al passo giusto, scambia individui
        if (enable_migration && size_ > 1 && gen > 0 && gen % migrate_every == 0) {
            migrationAllgather(200);
        }
    }

    // 2) Aggiorna best locale
    orderPopulation();
    best_path_    = population_[sorted_indices_[0]];
    best_fitness_ = fitness_[sorted_indices_[0]];

    // 3) Scrivi evolution.dat
    {
       string evo_file = "evolution_rank" + to_string(rank_) + suffix + ".dat";
        ofstream fevo(evo_file);
        fevo << "# Generation Best_Length AvgBestHalf\n";
        for (auto& d : evolution_data_) {
            fevo << d.generation << " "
                 << d.best_length << " "
                 << d.avg_best_half_length << "\n";
        }
        fevo.close();
    }

    // 4) Scrivi best_path.dat
    {
        string path_file = "best_path_rank" + to_string(rank_) + suffix + ".dat";
        ofstream fpath(path_file);
        for (int idx : best_path_) {
            const auto& c = cities_[idx];
            fpath << c.x << " " << c.y << "\n";
        }
        // Chiudi il percorso tornando alla prima città
        const auto& c0 = cities_[best_path_[0]];
        fpath << c0.x << " " << c0.y << "\n";
        fpath.close();
    }
}


// --- recupero best
double GA::getBestFitness() const {
    return best_fitness_;
}
// Getter per il percorso del migliore individuo locale
const vector<int>& GA::getBestPath() const {
    return best_path_;
}

// Implementazione della migrazione collettiva: scambia i migliori n_people individui
void GA::migrationAllgather(int n_people) {
    // 1. Prepara il buffer di invio con i migliori n_people individui
    vector<int> send_buffer(n_people * N_CITIES);
    auto top_indices = getTopIndices(n_people, true);
    
    for (int i = 0; i < n_people; ++i) {
        copy(population_[top_indices[i]].begin(),
                  population_[top_indices[i]].end(),
                  send_buffer.begin() + i * N_CITIES);
    }

    // 2. Esegui l'Allgather per scambiare gli individui tra tutti i nodi
    vector<int> all_buffer(size_ * n_people * N_CITIES);
    MPI_Allgather(
        send_buffer.data(),
        n_people * N_CITIES,
        MPI_INT,
        all_buffer.data(),
        n_people * N_CITIES,
        MPI_INT,
        MPI_COMM_WORLD
    );

    // 3. Seleziona casualmente n_people continenti diversi (escludendo il proprio rank)
    vector<int> selected_src;
    vector<int> available_src;
    for (int i = 0; i < size_; i++) {
        if (i != rank_) {
            available_src.push_back(i);
        }
    }
    
    // Mescola i rank disponibili usando rnd_
    for (int i = available_src.size() - 1; i > 0; --i) {
        int j = static_cast<int>(rnd_.Rannyu() * (i + 1));
        swap(available_src[i], available_src[j]);
    }
    
    // Prendi i primi min(n_people, available_src.size()) rank
    int num_immigrants = min(n_people, static_cast<int>(available_src.size()));
    if (num_immigrants > 0) {
        selected_src.assign(available_src.begin(), available_src.begin() + num_immigrants);
    }

    // 4. Estrai i migliori individui da ciascuna sorgente selezionata
    vector<vector<int>> immigrants;
    for (int src : selected_src) {
        int offset = src * n_people * N_CITIES;
        // Prendi il MIGLIORE individuo da questa sorgente
        vector<int> best_immigrant(N_CITIES);
        copy_n(all_buffer.begin() + offset, N_CITIES, best_immigrant.begin());
        immigrants.push_back(move(best_immigrant));
    }

    // 5. Sostituisci i peggiori individui nella popolazione corrente
    if (!immigrants.empty()) {
        auto worst_indices = getTopIndices(immigrants.size(), false);
        for (size_t i = 0; i < immigrants.size(); ++i) {
            population_[worst_indices[i]] = immigrants[i];
            fitness_[worst_indices[i]] = calculateFitness(immigrants[i]);
        }
    }
}

// --- Inizializzazione popolazione: genera POP_SIZE permutazioni valide
void GA::initializePopulation() {
    population_.resize(POP_SIZE, vector<int>(N_CITIES));
    fitness_.resize(POP_SIZE);
    
    for (int i = 0; i < POP_SIZE; ++i) {
        auto& p = population_[i];
        p[0] = 0;
        for (int j = 1; j < N_CITIES; ++j) p[j] = j;
        
        // Fisher-Yates shuffle (preserve starting city)
        for (int k = N_CITIES-1; k > 1; --k) {
            int j = 1 + int(rnd_.Rannyu() * k);
            swap(p[k], p[j]);
        }
        fitness_[i] = calculateFitness(p);
    }
}

// Un’evoluzione: selezione, crossover, mutazioni, elitismo
void GA::evolveOneGeneration() {
    orderPopulation();
    vector<vector<int>> new_population;
    vector<double> new_fitness;

    // Elitismo: conserva il milgiore
    new_population.push_back(population_[sorted_indices_[0]]);
    new_fitness.push_back(fitness_[sorted_indices_[0]]);

     // Genera nuovi individui fino a POP_SIZE
    while (new_population.size() < POP_SIZE) {
        // Selezione proporzionale (esponentiale)
        int p1 = sorted_indices_[min(POP_SIZE-1,
            static_cast<int>(POP_SIZE * pow(rnd_.Rannyu(), SELECTION_P)))];
        int p2 = sorted_indices_[min(POP_SIZE-1,
            static_cast<int>(POP_SIZE * pow(rnd_.Rannyu(), SELECTION_P)))];

        vector<int> child;
        if (rnd_.Rannyu() < CROSS_RATE) {
            child = crossoverOX(population_[p1], population_[p2]);
        } else {
            child = population_[rnd_.Rannyu() < 0.5 ? p1 : p2];
        }

        // applica mutazioni 
        if (rnd_.Rannyu() < MUT_RATE) swapMut(child);
        if (rnd_.Rannyu() < MUT_RATE) invertMut(child);
        if (rnd_.Rannyu() < MUT_RATE) shiftMut(child);
        if (rnd_.Rannyu() < MUT_RATE) edgeOptimizerMut(child);
        if (rnd_.Rannyu() < MUT_RATE) blockFlipMut(child);

        // Vector child è valido? Aggiungi al pool
        if (isValid(child)) {
            new_population.push_back(child);
            new_fitness.push_back(calculateFitness(child));
        }
    }
// Sostituisci popolazione e fitness
    population_ = move(new_population);
    fitness_    = move(new_fitness);
}




// Restituisce i primi n indici migliori (take_best=true) o peggiori (false)
vector<int> GA::getTopIndices(int n, bool take_best) {
    vector<int> indices(POP_SIZE);
    iota(indices.begin(), indices.end(), 0);
    
    sort(indices.begin(), indices.end(), [&](int a, int b) {
        return take_best ? (fitness_[a] > fitness_[b]) : (fitness_[a] < fitness_[b]);
    });
    
    return vector<int>(indices.begin(), indices.begin() + n);
}



// Ordina sorted_indices_ per fitness decrescente
void GA::orderPopulation() {
    sorted_indices_.resize(POP_SIZE);
    iota(sorted_indices_.begin(), sorted_indices_.end(), 0);
    sort(sorted_indices_.begin(), sorted_indices_.end(),
        [&](size_t i, size_t j) { return fitness_[i] > fitness_[j]; });
}
// Calcola fitness = 1 / lunghezza del percorso
double GA::calculateFitness(const vector<int>& p) const {
    double sum = 0.0;
    for (int i = 0; i < N_CITIES; ++i) {
        int j = (i+1) % N_CITIES;
        sum += dist(p[i], p[j]);
    }
    return 1.0 / sum;
}
// Mutazione locale: mischia un sottoblocco di 3–5 città
void GA::localScrambleMut(vector<int>& p) {
    int start = 1 + int(rnd_.Rannyu() * (N_CITIES - 3));
    int end   = start + 2 + int(rnd_.Rannyu() * 3);
    if (end >= N_CITIES) end = N_CITIES - 1;

    // Fisher–Yates shuffle solo fra gli indici [start .. end-1]
    for (int i = end - 1; i > start; --i) {
        int j = start + int(rnd_.Rannyu() * (i - start + 1));  
        swap(p[i], p[j]);
    }
}

// Verifica che p sia una permutazione valida che inizia da 0
bool GA::isValid(const vector<int>& p) const {
    if (p[0] != 0) return false;
    vector<bool> seen(N_CITIES, false);
    for (int city : p) {
        if (city < 0 || city >= N_CITIES || seen[city]) return false;
        seen[city] = true;
    }
    return true;
}
// Order Crossover (OX): prende un sottoarray da A e riempie il resto con B
vector<int> GA::crossoverOX(const vector<int>& A, const vector<int>& B) {
    int c1 = 1 + int(rnd_.Rannyu()*(N_CITIES-2));
    int c2 = 1 + int(rnd_.Rannyu()*(N_CITIES-2));
    if (c1 > c2) swap(c1, c2);

    vector<int> child(N_CITIES, -1);
    vector<bool> used(N_CITIES, false);
    
    // Copia il segmento [c1,c2] da A
    for(int i = c1; i <= c2; ++i) {
        child[i] = A[i];
        used[A[i]] = true;
    }

    // Riempimento da B dopo c2 (con wrap-around)
    int pos = (c2 + 1) % N_CITIES;  // Inizia subito dopo c2
    int start_B = pos;               // Punto di partenza in B
    
    for(int k = 0; k < N_CITIES; ++k) {
    int city = B[(start_B + k) % N_CITIES];
    if (!used[city]) {
        // Cerca la prossima posizione libera in ordine sequenziale
        while (child[pos % N_CITIES] != -1) {
            pos++;
        }
        child[pos % N_CITIES] = city;
        used[city] = true;
        pos++;
    }
}

    return child;
}
// Scambia due città a caso (swap mutation)
void GA::swapMut(vector<int>& p) {
    int a = 1 + int(rnd_.Rannyu()*(N_CITIES-1));
    int b = 1 + int(rnd_.Rannyu()*(N_CITIES-1));
    swap(p[a], p[b]);
}
// Inverte un sottoblocco (inversion mutation)
void GA::invertMut(vector<int>& p) {
    int a = 1 + int(rnd_.Rannyu()*(N_CITIES-1));
    int b = 1 + int(rnd_.Rannyu()*(N_CITIES-1));
    if (a>b) swap(a,b);
    reverse(p.begin()+a, p.begin()+b+1);
}

// Pre-calcola tutte le distanze Euclidee fra le città in una matrice
void GA::precomputeDistances() {
    distance_matrix_.resize(N_CITIES, vector<double>(N_CITIES));
    for (int i = 0; i < N_CITIES; ++i) {
        for (int j = 0; j < N_CITIES; ++j) {
            double dx = cities_[i].x - cities_[j].x;
            double dy = cities_[i].y - cities_[j].y;
            distance_matrix_[i][j] = sqrt(dx*dx + dy*dy);
        }
    }
}

// Ritorna la distanza fra due città usando la matrice
double GA::dist(int a, int b) const {
    return distance_matrix_[a][b];
}



// Sposta un blocco di lunghezza m di città più avanti di shift posizioni
void GA::shiftMut(vector<int>& p) {
    int m = 1 + int(rnd_.Rannyu() * (N_CITIES - 2)); // m <= N-2

    int max_start = N_CITIES - m - 1;
    if (max_start < 1) return;
    int start = 1 + int(rnd_.Rannyu() * (max_start - 1));

    int max_shift = N_CITIES - start - m;
    if (max_shift <= 0) return;
    int shift = int(rnd_.Rannyu() * (max_shift));

    vector<int> block(p.begin() + start, p.begin() + start + m);
    p.erase(p.begin() + start, p.begin() + start + m);
    p.insert(p.begin() + start + shift, block.begin(), block.end());
}
// Scambia due blocchi di lunghezza m
void GA::blockSwapMut(vector<int>& p) {
    if (N_CITIES < 4) return;

    // Dimensione massima del blocco
    int max_block_size = (N_CITIES - 1) / 2;
    if (max_block_size < 1) return;
    int m = 1 + int(rnd_.Rannyu() * (max_block_size - 1));

    // Calcolo di start1 e start2 con controlli di validità
    int max_start1 = N_CITIES - 2 * m - 1;
    if (max_start1 < 1) return;
    int start1 = 1 + int(rnd_.Rannyu() * max_start1);

    int remaining_space = N_CITIES - start1 - 2 * m;
    if (remaining_space <= 0) return;
    int start2 = start1 + m + int(rnd_.Rannyu() * remaining_space);

    // Swap dei blocchi
    for (int i = 0; i < m; ++i) {
        swap(p[start1 + i], p[start2 + i]);
    }
}
// Utility: carica da file le coordinate delle città
vector<City> loadCitiesFromFile(const string& filename) {
    vector<City> cities;
    ifstream fin(filename);
    double x, y;
    while (fin >> x >> y) {
        cities.push_back({x, y});
    }
    return cities;
}
// Ottimizzazione locale: scambia vicini se migliora
void GA::edgeOptimizerMut(vector<int>& p) {
    for (int i = 1; i < N_CITIES - 1; ++i) {
        double original = dist(p[i-1], p[i]) + dist(p[i], p[i+1]);
        swap(p[i], p[i+1]);
        double nuovo = dist(p[i-1], p[i]) + dist(p[i], p[i+1]);
        if (nuovo >= original) swap(p[i], p[i+1]); // Annulla se non migliora
    }
}

// Inversione di blocco di lunghezza variabile
void GA::blockFlipMut(vector<int>& p) {
    int start = 1 + int(rnd_.Rannyu() * (N_CITIES - 5));
    int end = start + 3 + int(rnd_.Rannyu() * 4);
    if (end >= N_CITIES) end = N_CITIES - 1;
    reverse(p.begin() + start, p.begin() + end + 1);
}