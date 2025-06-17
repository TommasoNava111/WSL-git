/*
 * Variational Monte Carlo per il ground state di una particella 1D nel potenziale V(x)=x^4 - (5/2)x^2.
 * 1) Definizione del trial wavefunction Ψ_T(x) = e^{−(x−μ)^2/(2σ^2)} + e^{−(x+μ)^2/(2σ^2)}, parametrizzato da σ e μ.
 * 2) Campionamento di |Ψ_T(x)|^2 con algoritmo di Metropolis: proposta uniforme di ampiezza δ.
 * 3) Calcolo dell’energia locale E_L(x) = −½ (Ψ_T''(x)/Ψ_T(x)) + V(x) usando la derivata seconda analitica.
 * 4) Metodo del data blocking: suddivisione dei M lanci in N blocchi per stimare media ed errore statistico.
 * 5) Stima finale di ⟨E⟩ e dell’acceptance ratio, ripetibile variando σ e μ per minimizzare l’energia.
 */
#include <iostream>   // IO basico (cout, cin, cerr)
#include <fstream>    // Gestione file (ifstream, ofstream)
#include <string>     // Stringhe in C++
#include <cmath>      // Funzioni matematiche (exp, pow, sqrt)
#include "random.h"   // Classe Random per generatore di numeri casuali

using namespace std;

// Legge due numeri primi da file per inizializzare il generatore
void ReadPrimes(const string& filename, int& p1, int& p2) {
    ifstream filePrimes(filename);
    if (!filePrimes.is_open()) { // Controllo apertura
        cerr << "ERRORE: Impossibile aprire il file " << filename << endl;
        exit(1);
    }
    filePrimes >> p1 >> p2;      // Carica i due primi
    filePrimes.close();
}

// Estrae semi da file e inizializza il generatore Random
void InitializeRandom(Random& rnd, const string& seedFile, int p1, int p2) {
    ifstream fileSeed(seedFile);
    if (!fileSeed.is_open()) {
        cerr << "ERRORE: Impossibile aprire il file " << seedFile << endl;
        exit(1);
    }
    string keyword;
    int seedArray[4];
    // Scorre il file alla ricerca di "RANDOMSEED"
    while (fileSeed >> keyword) {
        if (keyword == "RANDOMSEED") {
            // Legge quattro interi da usare come semi
            fileSeed >> seedArray[0] >> seedArray[1] >> seedArray[2] >> seedArray[3];
            rnd.SetRandom(seedArray, p1, p2);
            break;
        }
    }
    fileSeed.close();
}

// Ritorna |\Psi_T(x)|^2 per il trial wavefunction gaussiano bimodale
double psiT_squared(double x, double sigma, double mu) {
    // Due gaussiane centrate in +mu e -mu
    double term1 = exp(-pow(x - mu, 2) / (2 * sigma * sigma));
    double term2 = exp(-pow(x + mu, 2) / (2 * sigma * sigma));
    // Il quadrato della somma
    return pow(term1 + term2, 2);
}

// Calcola l'energia locale E_L(x) = (-1/2 · Ψ''/Ψ) + V(x)
double compute_energy_local(double x, double sigma, double mu) {
    // Ricostruisco term1 e term2 come sopra
    double term1 = exp(-pow(x - mu, 2) / (2 * sigma * sigma));
    double term2 = exp(-pow(x + mu, 2) / (2 * sigma * sigma));
    
    // Derivata seconda analitica di Ψ_T divisa per Ψ_T:
    // ( ( (x-μ)^2/σ^4 - 1/σ^2 ) e^{...} + idem con +μ ) / (e^{...}+e^{...})
    double numerator = ((pow(x - mu, 2)/(pow(sigma,4)) - 1.0/(sigma*sigma)) * term1)
                     + ((pow(x + mu, 2)/(pow(sigma,4)) - 1.0/(sigma*sigma)) * term2);
    double denominator = term1 + term2;
    double kinetic = -0.5 * (numerator / denominator);
    
    // Potenziale V(x) = x^4 - (5/2) x^2
    double potential = pow(x,4) - (5.0/2.0)*pow(x,2);
    
    return kinetic + potential;
}

// Esegue il ciclo di equilibrazione e poi il data blocking su M tossi
void ComputeEnergy(Random& rnd, double sigma, double mu, double delta,
                   int total_throws, int nBlocks, int nEquilibration,
                   ofstream& out, double& avg_energy, double& error_energy, double& acceptance) {
    int block_size = total_throws / nBlocks; // Tossic per blocco
    double x = 0.0;                          // Punto di partenza in 0
    
    // ===== Fase di equilibrazione Metropolis =====
    for (int i = 0; i < nEquilibration; ++i) {
        double x_new = x + delta * (2 * rnd.Rannyu() - 1.0); // Proposta uniforme in [-δ,+δ]
        double A = psiT_squared(x_new, sigma, mu) / psiT_squared(x, sigma, mu);
        if (rnd.Rannyu() < A) { // Accetta secondo criterio di Metropolis
            x = x_new;
        }
    }
    
    // ===== Campionamento a blocchi =====
    double cum_avg = 0.0, cum_avg2 = 0.0; // Sommatorie per errori
    int n_accepted = 0, n_trial = 0;      // Contatori accettazioni
    
    for (int block = 0; block < nBlocks; ++block) {
        double block_sum = 0.0;
        for (int j = 0; j < block_size; ++j) {
            // Proposta Metropolis
            double x_new = x + delta * (2 * rnd.Rannyu() - 1.0);
            double A = psiT_squared(x_new, sigma, mu) / psiT_squared(x, sigma, mu);
            n_trial++;
            if (rnd.Rannyu() < A) {
                x = x_new;
                n_accepted++;
            }
            // Accumula energia locale
            block_sum += compute_energy_local(x, sigma, mu);
        }
        // Media del blocco
        double block_avg = block_sum / block_size;
        cum_avg  += block_avg;
        cum_avg2 += block_avg * block_avg;

        // Calcola media cumulativa ed errore
        double avg   = cum_avg  / (block + 1);
        double avg2  = cum_avg2 / (block + 1);
        double error = (block > 0) ? sqrt((avg2 - avg*avg) / block) : 0.0;
        
        // Scrive su file: numero di tossi, media cumulativa, errore, acceptance parziale
        out << (block + 1) * block_size << " " << avg << " " << error << " "
            << static_cast<double>(n_accepted) / n_trial << endl;
    }

    // Ritorna valori globali
    avg_energy   = cum_avg  / nBlocks;
    error_energy = (nBlocks > 1) 
        ? sqrt((cum_avg2/nBlocks - avg_energy*avg_energy) / (nBlocks - 1)) 
        : 0.0;
    acceptance   = (n_trial > 0) ? static_cast<double>(n_accepted) / n_trial : 0.0;
}

int main(int argc, char *argv[]) {
    Random rnd;              // Istanza generatore RNG
    int p1, p2;              
    ReadPrimes("Primes", p1, p2);                      // Legge i primi
    InitializeRandom(rnd, "seed.in", p1, p2);           // Inizializza RNG
    
    // ----- Parametri variationali iniziali -----
    double sigma = 0.59;    // Larghezza gaussiana                      0.6          0.61
    double mu    = 0.82;    // Distanza dei due picchi dal centro       0.8    oppure 0.81 con 2.5 e 0.6 -0.455911 ± 0.00746134
    double delta = 2.5;    // Ampiezza proposta Metropolis               2.5        2.6
    int M     = 10000;     // Numero totale di lanci
    int N     = 100;       // Numero di blocchi
    int equilibration = 5000; // Tossic di equilibrazione
    
    // Prepara file di output
    ofstream out_energy("energy.txt");
    if (!out_energy) {
        cerr << "ERRORE: Impossibile aprire energy.txt" << endl;
        return 1;
    }
    out_energy << "throws avg error acceptance" << endl;
    
    double energy, error, acc;
    // Esegue il VMC e riempie energy.txt
    ComputeEnergy(rnd, sigma, mu, delta, M, N, equilibration, out_energy, energy, error, acc);
    
    // Stampa a video i risultati finali
    cout << "Energia media = "    << energy << " ± " << error << endl;
    cout << "Acceptance ratio globale = " << acc   << endl;
    
    rnd.SaveSeed();        // Salva lo stato del RNG
    out_energy.close();    // Chiude il file
    return 0;
}
