/*
  // 1) Inizializzazione
  //    - Partiamo da σ = 1.0 e μ = 2.0 (scelta arbitraria).
  //    - Impostiamo temperatura iniziale T = T0 (es. 6.0).
  //    - Definiamo il passo di proposta per i parametri: Δparam = param_step (es. 0.05).
  //
  // 2) Loop di Simulated Annealing (nsteps iterazioni):
  //    a) Valutazione corrente:
  //       • Chiamo vmc(rnd, σ, μ, vmc_steps, delta) per ottenere:
  //           – energy_curr  = energia media stimata per (σ, μ).
  //           – error_curr   = incertezza statistica.
  //           – accept_rate  = tasso di accettanza del Metropolis interno.
  //
  //    b) Generazione di una proposta (σ_new, μ_new):
  //       • σ_new = σ + (2*U(0,1) - 1) * Δparam
  //       • μ_new = μ + (2*U(0,1) - 1) * Δparam
  //       • Se σ_new ≤ 0.1 o μ_new ≤ 0.1, scartiamo la proposta e manteniamo i valori vecchi.
  //
  //    c) Valutazione proposta:
  //       • Chiamo vmc(rnd, σ_new, μ_new, vmc_steps, delta) per ottenere energy_prop.
  //
  //    d) Criterio di accettazione Metropolis–SA:
  //       • Calcolo ΔE = energy_prop - energy_curr.
  //       • Se ΔE < 0 ⇒ accetto sempre la proposta (σ = σ_new, μ = μ_new).
  //       • Altrimenti ⇒ accetto con probabilità p = exp(-ΔE / T).
  //           – Genero r = U(0,1); se r < p accetto, altrimenti riffuto.
  //
  //    e) Registrazione:
  //       • Salvo su file “sa_data.txt”: passo, energy_curr, error_curr, σ, μ, accept_rate.
  //
  //    f) Aggiornamento temperatura:
  //       • T = T * cool   (es. cool = 0.95)
  //
  // 3) Fine ciclo SA:
  //    - Il file “sa_data.txt” contiene la storia di (σ, μ) e delle energie.
  //
  // 4) Estrazione dei parametri ottimali:
  //    - find_optimal_parameters legge “sa_data.txt” riga per riga.
  //    - Tiene traccia del valore minimo di energy_curr.
  //    - Restituisce (σ_best, μ_best) corrispondenti a quell’energia minima.
  //
  // 5) Misura finale con parametri ottimali:
  //    - Chiamo vmc(rnd, σ_best, μ_best, nsteps_final, delta_final, save_blocks=true, save_samples=true)
  //      per ottenere stima definitiva di energia, errore e campioni x.
  //    - Salvo risultati in “final_energy.txt” e campioni in “psi_histogram.txt”.
  //
  // In questo modo l’algoritmo esplora sistematicamente lo spazio (σ, μ),
  // bilanciando esplorazione (accettanze casuali) e sfruttamento (accetto sempre ΔE<0),
  // e riduce gradualmente la temperatura per convergere verso il minimo globale.
*/


#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <limits>
#include "random.h"

using namespace std;

// Physical constants: reduced Planck's constant and mass set to 1 for convenience
const double hbar = 1.0;
const double m = 1.0;

// ##################################################
// ## Funzioni di inizializzazione RNG
// ##################################################

// Legge due numeri primi da file per inizializzare RNG
void ReadPrimes(const string& filename, int& p1, int& p2) {
    ifstream filePrimes(filename);                      // Apre file in lettura
    if (!filePrimes) {                                  // Controlla apertura
        cerr << "Errore apertura " << filename << endl;
        exit(1);                                        // Termina in caso di errore
    }
    filePrimes >> p1 >> p2;                             // Estrae due interi primi
    filePrimes.close();                                 // Chiude il file
}

// Inizializza RNG leggendo seme da file seedFile e numeri primi p1,p2
void InitializeRandom(Random& rnd, const string& seedFile, int p1, int p2) {
    ifstream fileSeed(seedFile);
    string keyword;
    int seed[4];                                        // Array per 4 semi interi
    
    // Scorre file fino a trovare la parola "RANDOMSEED"
    while(fileSeed >> keyword) {
        if(keyword == "RANDOMSEED") {
            fileSeed >> seed[0] >> seed[1] >> seed[2] >> seed[3];  // Legge i semi
            rnd.SetRandom(seed, p1, p2);                         // Imposta RNG con semi e primi
            break;                                                // Esce dal loop
        }
    }
    fileSeed.close();                                // Chiude file seme
}

// ##################################################
// ## Funzioni VMC
// ##################################################

// Valuta il trial wave function Ψ_T(x;σ,μ)
double psi(double x, double sigma, double mu) {
    // Somma di due gaussiane centrata in +μ e -μ con varianza σ²
    return exp(-pow(x-mu,2)/(2*sigma*sigma)) + exp(-pow(x+mu,2)/(2*sigma*sigma));
}

// Calcolo della derivata seconda logaritmica di Ψ_T, usata nell'energia locale
double seconda_derivata(double x, double sigma, double mu) {
    // Termini esponenziali separati (gaussiana dx)
    double t1 = exp(-pow(x-mu,2)/(2*sigma*sigma));
    double t2 = exp(-pow(x+mu,2)/(2*sigma*sigma));
    
    // Contributi della derivata seconda di ogni gaussiana: (x-μ)²/σ^4 - 1/σ^2
    double d1 = (pow(x-mu,2)/pow(sigma,4) - 1/pow(sigma,2)) * t1;
    double d2 = (pow(x+mu,2)/pow(sigma,4) - 1/pow(sigma,2)) * t2;
    
    // Somma e normalizzazione dividendo per Ψ(x)
    return (d1 + d2)/(t1 + t2);
}

// Calcola l'energia locale E_L(x) = -½(Ψ''/Ψ) + V(x)
double energia_locale(double x, double sigma, double mu) {
    // Parte cinetica con hbar^2/(2m)=½ e potenziale V(x)=x^4-2.5 x^2
    return -0.5*seconda_derivata(x, sigma, mu) + pow(x,4) - 2.5*x*x;
}

// Struttura per restituire i risultati di una singola run VMC
struct VMCResult {
    double energy;            // Energia media stimata
    double error;             // Errore statistico (std. error)
    double accept_rate;       // Frazione di passi accettati
    vector<double> samples;   // Campioni x salvati (opzionale)
    double final_delta;       // Passo delta finale dopo termalizzazione
};

// Implementazione del algoritmo VMC con Metropolis
VMCResult vmc(Random& rnd,
              double sigma,
              double mu,
              int nsteps,
              double delta,
              bool save_blocks=false,
              bool save_samples=false) {

    vector<double> energies;               // Energia locale ad ogni passo
    vector<double> samples;
    energies.reserve(nsteps);
    if (save_samples) samples.reserve(nsteps);

    double x = 0.0;                        // Posizione iniziale
    int accepted_therm = 0;                // Contatore accettati in termalizzazione
    int accepted_meas  = 0;                // Contatore accettati in misurazione

    // --- Termalizzazione con adattamento di delta per accettanza ~50% ---
    const int n_therm = 5000;
    const int adjust_interval = 100;      // Intervallo per aggiustare delta
    double target_accept = 0.5;
    int accepted_in_interval = 0;

    for (int i = 0; i < n_therm; ++i) {
        double x_new = x + (2*rnd.Rannyu() - 1) * delta;                // Proposta uniforme
        double A = pow(psi(x_new, sigma, mu) / psi(x, sigma, mu), 2);   // Rapporto di probabilità
        if (rnd.Rannyu() < A) {
            x = x_new;
            accepted_therm++;
            accepted_in_interval++;
        }
        // Ogni adjust_interval passi, aggiusta delta
        if ((i + 1) % adjust_interval == 0) {
            double current_accept = double(accepted_in_interval) / adjust_interval;
            if (current_accept > target_accept)
                delta *= 1.1;
            else
                delta *= 0.9;
            // Limita delta in [0.1,5.0]
            delta = max(0.1, min(delta, 5.0));
            accepted_in_interval = 0;
        }
    }

    // --- Fase di misurazione vera e propria ---
    for (int i = 0; i < nsteps; ++i) {
        double x_new = x + (2*rnd.Rannyu() - 1) * delta;
        double A = pow(psi(x_new, sigma, mu) / psi(x, sigma, mu), 2);
        if (rnd.Rannyu() < A) {
            x = x_new;
            accepted_meas++;
        }
        double e = energia_locale(x, sigma, mu);  // Calcola energia locale
        energies.push_back(e);
        if (save_samples) samples.push_back(x);
    }

    // --- Data blocking per stima dell'errore ---
    const int nblocks = 100;
    int L = nsteps / nblocks;
    vector<double> block_averages;
    block_averages.reserve(nblocks);
    for(int b=0; b<nblocks; ++b) {
        double block_sum = 0.0;
        for(int j=0; j<L; ++j) {
            block_sum += energies[b*L + j];
        }
        block_averages.push_back(block_sum / L);
    }
    // Calcolo media e errore su blocchi
    double sum = 0.0, sum2 = 0.0;
    for(double avg : block_averages) {
        sum += avg;
        sum2 += avg*avg;
    }
    double energy = sum / nblocks;
    double error = sqrt((sum2/nblocks - energy*energy) / (nblocks - 1));

    // Salvataggio opzionale dei blocchi su file
    if(save_blocks) {
        ofstream block_file("final_energy_blocks.txt");
        double running_sum = 0.0, running_sum2 = 0.0;
        for(int b=0; b<nblocks; ++b) {
            running_sum += block_averages[b];
            running_sum2 += block_averages[b]*block_averages[b];
            double mean = running_sum / (b+1);
            double err = sqrt((running_sum2/(b+1) - mean*mean) / (b+1));
            block_file << (b+1) << " " << mean << " " << err << endl;
        }
    }

    // Ritorna strutt VMCResult con risultati e campioni
    return { energy, error, double(accepted_meas)/nsteps, samples, delta };
}

// ##################################################
// ## Simulated Annealing
// ##################################################

// Struttura per leggere i dati SA e trovare la miglior configurazione
struct SAState {
    int step;             // Passo SA
    double energy;        // Energia corrispondente
    double sigma;         // Parametro σ
    double mu;            // Parametro μ
    double accept_rate;   // Accettanza in VMC
};

// Cerca nel file SA il minimo di energia e restituisce parametri ottimi
SAState find_optimal_parameters(const string& filename) {
    ifstream sa_file(filename);
    SAState optimal = {0, numeric_limits<double>::max(), 0, 0, 0};
    int step;
    double energy, err, sigma, mu, accept_rate;
    while(sa_file >> step >> energy >> err >> sigma >> mu >> accept_rate) {
        if(energy < optimal.energy) {
            optimal = {step, energy, sigma, mu, accept_rate};
        }
    }
    return optimal;
}

// Simulated Annealing per ottimizzare (σ,μ)
void simulated_annealing(Random& rnd,
                        int nsteps,
                        double T0,
                        double cool,
                        double param_step,
                        int vmc_steps,
                        double delta) {
    ofstream sa_file("sa_data.txt");       // File output dati SA
    double sigma = 1.0, mu = 2.0;           // Condizione iniziale
    double T = T0;                          // Temperatura iniziale

    for(int step=0; step<nsteps; step++) {
        // Calcola energia corrente nel punto (σ,μ)
        VMCResult current = vmc(rnd, sigma, mu, vmc_steps, delta);
        
        // Campiona nuova proposta di (σ,μ)
        double sigma_new = sigma + (2*rnd.Rannyu()-1)*param_step;
        double mu_new    = mu    + (2*rnd.Rannyu()-1)*param_step;
        // Mantiene σ,μ >0.1
        if(sigma_new <= 0.1 || mu_new <= 0.1) { 
            sigma_new = sigma;
            mu_new    = mu;
        }

        // Valuta proposta
        VMCResult proposal = vmc(rnd, sigma_new, mu_new, vmc_steps, delta);

        double dE = proposal.energy - current.energy;
        // Accetto se migliora o con probabilità exp(-dE/T)
        if(dE < 0 || exp(-dE/T) > rnd.Rannyu()) {
            sigma = sigma_new;
            mu    = mu_new;
        }

        // Registra dati: passo, energia, errore, σ, μ, accettanza
        sa_file << step << " "
                << current.energy << " "
                << current.error  << " "
                << sigma          << " "
                << mu             << " "
                << current.accept_rate << endl;
        
        // Raffredda la temperatura
        T *= cool;
    }
}

// ##################################################
// ## Main program
// ##################################################
int main() {
    Random rnd;
    int p1, p2;
    ReadPrimes("Primes", p1, p2);                           // Inizializza RNG
    InitializeRandom(rnd, "seed.in", p1, p2);

    // Esegue annealing con 200 passi, T0=6.0, fattore raffreddamento=0.95
    simulated_annealing(rnd, 200, 6.0, 0.95, 0.05, 10000, 2.0);

    // Trova parametri ottimali da file
    SAState optimal = find_optimal_parameters("sa_data.txt");
    cout << "Parametri ottimali:\n"
         << "σ = " << optimal.sigma << "\nμ = " << optimal.mu
         << "\nE = " << optimal.energy
         << "\nAccettanza: " << optimal.accept_rate << endl;

    // Misura finale con milioni di passi, salva blocchi e campioni
    VMCResult final = vmc(rnd, optimal.sigma, optimal.mu, 1000000, 1.0, true, true);
    
    // Salva risultato finale
    ofstream final_file("final_energy.txt");
    final_file << optimal.sigma << " " << optimal.mu << " "
               << final.energy << " " << final.error << " "
               << final.accept_rate << endl;

    // Salva ogni campione x in file per istogramma
    ofstream hist_file("psi_histogram.txt");
    for(double x : final.samples)
        hist_file << x << endl;

    rnd.SaveSeed();  // Salva seme per riproducibilità
    return 0;
}
