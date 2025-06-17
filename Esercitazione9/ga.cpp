/*
  ga.cpp

  Implementazione della classe GA per la soluzione del Traveling Salesman Problem
  mediante un algoritmo genetico. Contiene:

  - Costruttore GA(Random& rnd, bool on_circle): inizializza il generatore di numeri casuali,
    la modalità di posizionamento delle città (su circonferenza o nel quadrato) e il parametro
    di selezione.

  - run(): coordina l’intero processo evolutivo, genera le città, inizializza la popolazione,
    esegue le generazioni (selezione, crossover, mutazioni, elitismo), salva i risultati
    (lunghezza del percorso migliore e media della top-half) e il percorso ottimale finale.

  - orderPopulation(): ordina gli indici della popolazione in base alla fitness decrescente.

  - generateCities(): genera le coordinate delle N_CITIES su cerchio unitario o in un quadrato
    di area π.

  - initializePopulation(): crea una popolazione di permutazioni valide con Fisher–Yates,
    mantenendo la città 0 come punto di partenza.

  - evolveOneGeneration(): applica elitismo, selezione (metodo a potenza), crossover Order-X
    e una serie di mutazioni (swap, inversione, shift, block-swap, local scramble), controlla
    la validità dei figli e aggiorna la popolazione.

  - calculateFitness(): calcola la fitness come 1/(lunghezza del percorso) usando L2 o L1.

  - isValid(): verifica che ogni percorso contenga tutte le città esattamente una volta e
    che la città 0 sia all’inizio.

  - crossoverOX(): implementa un crossover che preserva l’ordine e i vincoli della permutazione.

  - Vari operatori di mutazione:
      * swapMut(): scambio di due città
      * invertMut(): inversione di un segmento
      * shiftMut(): spostamento di un blocco contiguo
      * blockSwapMut(): scambio di due blocchi
      * localScrambleMut(): mescolamento locale di un sottoinsieme

  - dist(): calcola la distanza tra due città (L2 o L1 a seconda della costante USE_L2).

  Insieme, queste funzioni realizzano un GA completo per ottimizzare il percorso del TSP
  su 34 città.
*/

#include "ga.h"  // Include la dichiarazione della classe GA
#include <algorithm> // Per sort, shuffle, iota ecc.
#include <numeric>  // Per iota
#include <fstream> // Per operazioni su file
#include <cmath>  // Per funzioni matematiche
using namespace std;
// Costruttore: riceve il RNG e un flag per generazione città su cerchio o quadrato
GA::GA(Random& rnd, bool on_circle) 
  : rnd_(rnd),  // Inizializza il riferimento al generatore
  on_circle_(on_circle), // Imposta il tipo di distribuzione delle città
  selection_p_(4.0) // Esponente di selezione predefinito
{
    suffix_ = on_circle_ ? "_circle" : "_square"; // Sceglie il suffisso per i file di output
}
// Ordina la popolazione in base alla fitness (decrescente)
void GA::orderPopulation() {
    // Crea vettore di indici
    sorted_indices_.resize(POP_SIZE); // Riserva spazio
    iota(sorted_indices_.begin(), sorted_indices_.end(), 0);  // Riempie con 0,1,...,POP_SIZE-1
    
     // Ordina gli indici in base al valore di fitness_, dal più grande al più piccolo
    sort(sorted_indices_.begin(), sorted_indices_.end(),
        [this](size_t i, size_t j) { return fitness_[i] > fitness_[j]; });
}
// Funzione principale che esegue il GA e salva i risultati
void GA::run() {
    generateCities();  // Genera le coordinate delle città
    initializePopulation();  // Inizializza la popolazione iniziale

    ofstream fout("results" + suffix_ + ".txt"); // File per best e avg per generazione
    ofstream fcity("cities" + suffix_ + ".txt"); // File con coordinate città
    ofstream fbest("best_path" + suffix_ + ".txt"); // File con percorso migliore finale

    // Salva le coordinate delle città
    for (auto& c : cities_) {
        fcity << c.x << " " << c.y << "\n";
    }
    fcity.close();

    // Ciclo di evoluzione
    for (int gen = 0; gen < GENERATIONS; ++gen) {
        evolveOneGeneration();  // Selezione → crossover → mutazioni

        // Copia e ordina la fitness in ordine decrescente
        vector<double> sorted_fitness = fitness_;
        sort(sorted_fitness.rbegin(), sorted_fitness.rend()); // dal migliore (valore più alto) al peggiore

        // Lunghezza del percorso migliore (inverto la fitness)
        double best_length = 1.0 / sorted_fitness[0];

        // Calcola la media delle lunghezze dei migliori 250 percorsi
        int top_n = 250; // metà della popolazione
        double avg_top_length = 0.0;
        for (int i = 0; i < top_n; ++i) {
            avg_top_length += 1.0 / sorted_fitness[i];
        }
        avg_top_length /= top_n;
 // Salva gen, best_length e avg_top_length su file
        fout << gen << " " << best_length << " " << avg_top_length << "\n";
    }
    fout.close();

    // Trova l'indice del miglior individuo (usando fitness_ non ordinata)
    int ib = max_element(fitness_.begin(), fitness_.end()) - fitness_.begin();
    vector<int> best_path = population_[ib];

    // Salva il percorso ottimizzato, ripetendo la città 0 alla fine
    for (int city : best_path) {
        fbest << city << "\n";
    }
    fbest << best_path[0] << "\n";
    fbest.close();
}
// Mutazione “local scramble”: mescola un sottoinsieme contiguo di città
void GA::localScrambleMut(vector<int>& p) {
    int start = 1 + int(rnd_.Rannyu() * (N_CITIES - 3));  // Inizio blocco
    int end   = start + 2 + int(rnd_.Rannyu() * 3);  // Fine blocco
    if (end >= N_CITIES) end = N_CITIES - 1; // Limite max

    // Fisher–Yates shuffle solo fra gli indici [start .. end-1]
    for (int i = end - 1; i > start; --i) {
        int j = start + int(rnd_.Rannyu() * (i - start + 1));  
        swap(p[i], p[j]);
    }
}
// Genera le coordinate delle città
void GA::generateCities() {
    cities_.resize(N_CITIES);
    if (on_circle_) {
        double step = 2*M_PI/N_CITIES; // Angolo tra le città
        for (int i = 0; i < N_CITIES; ++i) {
            double a = i*step;
            cities_[i] = { cos(a), sin(a) };  // Coordinate sul cerchio unitario
        }
    } else {
        double L = sqrt(M_PI); // Lato del quadrato con area π (come la circonferenza)
        for (int i = 0; i < N_CITIES; ++i)
            cities_[i] = { rnd_.Rannyu()*L, rnd_.Rannyu()*L };  // Coordinate uniformi nel quadrato
    }
}
// Inizializza la popolazione con permutazioni casuali valide
void GA::initializePopulation() {
    population_.assign(POP_SIZE, vector<int>(N_CITIES));
    fitness_.resize(POP_SIZE);
    for (int i = 0; i < POP_SIZE; ++i) {
        auto& p = population_[i];
        p[0] = 0;  // Fissa la città iniziale
        for (int j = 1; j < N_CITIES; ++j) p[j] = j; // Sequenza ordinata
        // Fisher-Yates (preserva indice 0)
        for (int k = N_CITIES-1; k > 1; --k) {
            int j = 1 + int(rnd_.Rannyu()*k);
            swap(p[k], p[j]);
        }
        // Calcola la fitness dell’individuo
        fitness_[i] = calculateFitness(p);
    }
}
// Costruisce la nuova generazione (elitismo + selezione + crossover + mutazioni)
void GA::evolveOneGeneration() {
    orderPopulation();  // Ordina la popolazione
    
    vector<vector<int>> new_population;
    vector<double> new_fitness;
    
    // Elitismo: mantieni il miglior individuo
    new_population.push_back(population_[sorted_indices_[0]]);
    new_fitness.push_back(fitness_[sorted_indices_[0]]);
 // Riempie il resto della popolazione
    while (new_population.size() < POP_SIZE) {
        // Selezione genitore 1
        double r1 = rnd_.Rannyu();
        int j1 = static_cast<int>(POP_SIZE * pow(r1, selection_p_));
        j1 = max(0, min(j1, POP_SIZE-1));
        int p1 = sorted_indices_[j1];

        // Selezione genitore 2
        double r2 = rnd_.Rannyu();
        int j2 = static_cast<int>(POP_SIZE * pow(r2, selection_p_));
        j2 = max(0, min(j2, POP_SIZE-1));
        int p2 = sorted_indices_[j2];

        // Crossover e mutazione (mantieni la logica esistente)
        vector<int> child;
        if(rnd_.Rannyu() < CROSS_RATE) {
            child = crossoverOX(population_[p1], population_[p2]);
        } else {
            // Altrimenti copia a caso p1 o p2
            child = (rnd_.Rannyu() < 0.5) ? population_[p1] : population_[p2];
        }

        // Applica mutazioni
        if(rnd_.Rannyu() < MUT_RATE) swapMut(child);
        if(rnd_.Rannyu() < MUT_RATE) invertMut(child);
        if(rnd_.Rannyu() < MUT_RATE) shiftMut(child);
        if(rnd_.Rannyu() < MUT_RATE) blockSwapMut(child);
        if(rnd_.Rannyu() < MUT_RATE) localScrambleMut(child);
// Verifica validità e inserisce in nuova popolazione
        if(isValid(child)) {
            new_population.push_back(child);
            new_fitness.push_back(calculateFitness(child));
        }
    }
// Sostituisce popolazione e fitness correnti
    population_.swap(new_population);
    fitness_.swap(new_fitness);
}
// Calcola la fitness come 1 / lunghezza del percorso
double GA::calculateFitness(const vector<int>& p) const {
    double sum = 0.0;
    for (int i = 0; i < N_CITIES; ++i) {
        int j = (i+1) % N_CITIES; // Visita circolare
        sum += dist(p[i], p[j]);
    }
    return 1.0/sum;
}
// Controlla che il percorso sia valido (ogni città esattamente una volta, p[0]==0)
bool GA::isValid(const vector<int>& p) const {
    if (p[0] != 0) return false;  // Prima città deve essere 0
    vector<bool> seen(N_CITIES,false);
    for (int x : p) {
        if (x<0||x>=N_CITIES||seen[x]) return false;
        seen[x] = true;
    }
    return true;
}
// Crossover Order X-preserving: copia un blocco da A e riempie da B nell’ordine
vector<int> GA::crossoverOX(const vector<int>& A, const vector<int>& B) {
    int c1 = 1 + int(rnd_.Rannyu()*(N_CITIES-2));
    int c2 = 1 + int(rnd_.Rannyu()*(N_CITIES-2));
    if (c1>c2) swap(c1,c2);

    vector<int> child(N_CITIES,-1);
    vector<bool> used(N_CITIES,false);
    child[0] = 0; used[0] = true;
 // Copia il segmento [c1..c2] da A
    for (int i = c1; i <= c2; ++i) {
        child[i] = A[i];
        used[A[i]] = true;
    }
     // Riempi il resto con i valori di B in ordine
    int pos = (c2+1) % N_CITIES;
    for (int k = 0; k < N_CITIES; ++k) {
        int city = B[(c2+1+k) % N_CITIES];
        if (!used[city]) {
            child[pos] = city;
            used[city] = true;
            pos = (pos+1) % N_CITIES;
        }
    }
    return child;
}
// Scambia due città a caso (eccetto la prima)
void GA::swapMut(vector<int>& p) {
    int a = 1 + int(rnd_.Rannyu()*(N_CITIES-1));
    int b = 1 + int(rnd_.Rannyu()*(N_CITIES-1));
    swap(p[a], p[b]);
}
// Inverte l'ordine di un segmento casuale (eccetto la prima città)
void GA::invertMut(vector<int>& p) {
    int a = 1 + int(rnd_.Rannyu()*(N_CITIES-1));
    int b = 1 + int(rnd_.Rannyu()*(N_CITIES-1));
    if (a>b) swap(a,b);
    reverse(p.begin()+a, p.begin()+b+1);
}
// Calcola la distanza L2 o L1 tra due città
double GA::dist(int a, int b) const {
    double dx = cities_[a].x - cities_[b].x;
    double dy = cities_[a].y - cities_[b].y;
    return USE_L2 ? (dx*dx + dy*dy) : (fabs(dx) + fabs(dy)); // Potenza due o Distanza L1
}


// Mutazione “shift”: sposta un blocco contiguo di m città in avanti di shift posizioni
void GA::shiftMut(vector<int>& p) {
    int m = 1 + int(rnd_.Rannyu() * (N_CITIES - 2)); // Dimensione blocco < N-1
    int start = 1 + int(rnd_.Rannyu() * (N_CITIES - m - 1));
    int shift = 1 + int(rnd_.Rannyu() * (N_CITIES - start - m - 1));
    
    // Sposta il blocco di m città di "shift" posizioni
    vector<int> block(p.begin() + start, p.begin() + start + m);
    p.erase(p.begin() + start, p.begin() + start + m);
    p.insert(p.begin() + start + shift, block.begin(), block.end());
}
// Mutazione “block swap”: scambia due blocchi di dimensione m
void GA::blockSwapMut(vector<int>& p) {
    if (N_CITIES < 4) return; // Almeno 2 blocchi da 2 città
    int m = 1 + int(rnd_.Rannyu() * ((N_CITIES-1)/2 - 1)); // m < N/2
    int start1 = 1 + int(rnd_.Rannyu() * (N_CITIES - 2*m - 1));
    int start2 = start1 + m + int(rnd_.Rannyu() * (N_CITIES - start1 - 2*m));
    
     // Scambia elemento per elemento
    for (int i = 0; i < m; ++i) {
        swap(p[start1 + i], p[start2 + i]);
    }
}