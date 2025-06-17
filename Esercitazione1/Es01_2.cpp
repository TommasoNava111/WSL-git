// Es01.2.cpp
//
// Questo programma estende il generatore di numeri pseudo-casuali Rannyu() come richiesto
// dall'esercizio 01.2, verificando il Teorema del Limite Centrale per tre "dice":
// 1) Dice standard (uniforme discreta 1–6)
// 2) Dice esponenziale p(x)=λ e^{-λx}, inversione cumulativa (λ=1)
// 3) Dice di Cauchy-Lorentz p(x)=1/π · [Γ/((x-μ)^2+Γ^2)], inversione cumulativa (μ=0, Γ=1)
//
// Struttura principale:
// - ReadPrimes(seedFile, p1, p2): legge due primi per l'inizializzazione RNG.
// - InitializeRandom(rnd, seedFile, p1, p2): carica quattro semi da seed.in alla keyword "RANDOMSEED".
// - EXPdist(λ, r): trasforma r∈(0,1) in X esponenziale via X = -1/λ·ln(r).
// - Cauchydist(μ, Γ, r): trasforma r∈(0,1) in Y di Cauchy via Y = μ + Γ·tan[π(r−0.5)].
// - ComputeUniformAverages(rnd, M, N, filename): per ogni N∈{1,2,10,100} calcola M=10^4 medie
//   campionarie S_N=(1/N)Σ_{i=1}^N r_i e salva su file.
// - ComputeExpAverages(rnd, M, N, λ, filename): analogamente con variabili esponenziali.
// - ComputeCauchyAverages(rnd, M, N, μ, Γ, filename): analogamente con variabili di Cauchy.
//
// Variabili di simulazione in main():
//   M = 10^4: numero di realizzazioni per ciascun N.
//   Ns = {1,2,10,100}: valori di N (numero di addendi in S_N).
//   λ = 1.0, μ = 0.0, Γ = 1.0: parametri delle distribuzioni continue.
//   fStd, fExp, fCau: nomi dei file di output per ciascun tipo di dice e ogni N.
//
// Output:
//   Per ogni N e per ciascuno dei tre metodi, un file .txt contenente M valori di S_N,
//   pronti per essere rappresentati come istogrammi in tre figure distinte.

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;

// -----------------------------------------------------------------------------
// ReadPrimes: legge due numeri primi dal file indicato e li restituisce via riferimento
// Se il file non si apre, stampa errore e termina con exit(1)
// -----------------------------------------------------------------------------
void ReadPrimes(const string& filename, int& p1, int& p2) {
    ifstream filePrimes(filename);
    if (!filePrimes.is_open()) {
        cerr << "ERRORE: Impossibile aprire il file " << filename << endl;
        exit(1);
    }
    filePrimes >> p1 >> p2;  // leggiamo i due interi
    filePrimes.close();
}

// -----------------------------------------------------------------------------
// InitializeRandom: legge il seme dal file seedFile e inizializza il generatore
// con i primi p1, p2. Cerca la keyword "RANDOMSEED" e poi quattro interi di seme.
// Se il file non si apre, stampa errore e termina con exit(1)
// -----------------------------------------------------------------------------
void InitializeRandom(Random& rnd, const string& seedFile, int p1, int p2) {
    ifstream fileSeed(seedFile);
    if (!fileSeed.is_open()) {
        cerr << "ERRORE: Impossibile aprire il file " << seedFile << endl;
        exit(1);
    }
    string keyword;
    int seedArr[4];
    // Leggiamo riga per riga fino a trovare RANDOMSEED
    while (fileSeed >> keyword) {
        if (keyword == "RANDOMSEED") {
            fileSeed >> seedArr[0] >> seedArr[1] >> seedArr[2] >> seedArr[3];
            rnd.SetRandom(seedArr, p1, p2);  // inizializzazione RNG
            break;
        }
    }
    fileSeed.close();
}

// -----------------------------------------------------------------------------
// EXPdist: trasforma un numero uniforme r in una variabile esponenziale
// Parametri:
//   lambda - parametro di scala (tasso) della distribuzione
//   r      - numero uniforme in (0,1)
// Ritorna X = -1/lambda * ln(r)
// -----------------------------------------------------------------------------
double EXPdist(double lambda, double r) {
    return -1.0/lambda * log(r);
}

// -----------------------------------------------------------------------------
// Cauchydist: trasforma un numero uniforme r in una variabile di Cauchy
// Parametri:
//   mu    - parametro di posizione (centro)
//   gamma - parametro di scala
//   r     - numero uniforme in (0,1)
// Ritorna X = mu + gamma * tan(pi*(r-0.5))
// -----------------------------------------------------------------------------
double Cauchydist(double mu, double gamma, double r) {
    return mu + gamma * tan(M_PI * (r - 0.5));
}

// -----------------------------------------------------------------------------
// ComputeUniformAverages:
//   - Genera M realizzazioni di S_N = (1/N)*sum_{i=1}^N r_i
//   - r_i sono numeri uniformi [0,1)
//   - Scrive le medie su file di testo, una per riga
// -----------------------------------------------------------------------------
void ComputeUniformAverages(Random& rnd, int M, int N, const string& filename) {
    ofstream out(filename.c_str());
    if (!out.is_open()) {
        cerr << "ERRORE: Impossibile aprire il file " << filename << endl;
        exit(1);
    }
    for (int i = 0; i < M; ++i) {
        double sum = 0.0;
        for (int k = 0; k < N; ++k) {
            sum += rnd.Rannyu();  // somma dei numeri uniformi
        }
        out << (sum / N) << "\n";  // media campionaria
    }
    out.close();
}

// -----------------------------------------------------------------------------
// ComputeExpAverages:
//   - Genera M realizzazioni di S_N = (1/N)*sum_{i=1}^N X_i
//   - X_i sono esponenziali con parametro lambda
//   - Scrive le medie su file di testo, una per riga
// -----------------------------------------------------------------------------
void ComputeExpAverages(Random& rnd, int M, int N, double lambda, const string& filename) {
    ofstream out(filename.c_str());
    if (!out.is_open()) {
        cerr << "ERRORE: Impossibile aprire il file " << filename << endl;
        exit(1);
    }
    for (int i = 0; i < M; ++i) {
        double sum = 0.0;
        for (int k = 0; k < N; ++k) {
            double r = rnd.Rannyu();
            sum += EXPdist(lambda, r);
        }
        out << (sum / N) << "\n";
    }
    out.close();
}

// -----------------------------------------------------------------------------
// ComputeCauchyAverages:
//   - Genera M realizzazioni di S_N = (1/N)*sum_{i=1}^N Y_i
//   - Y_i sono di Cauchy con parametri mu e gamma
//   - Scrive le medie su file di testo, una per riga
// -----------------------------------------------------------------------------
void ComputeCauchyAverages(Random& rnd, int M, int N, double mu, double gamma, const string& filename) {
    ofstream out(filename.c_str());
    if (!out.is_open()) {
        cerr << "ERRORE: Impossibile aprire il file " << filename << endl;
        exit(1);
    }
    for (int i = 0; i < M; ++i) {
        double sum = 0.0;
        for (int k = 0; k < N; ++k) {
            double r = rnd.Rannyu();
            sum += Cauchydist(mu, gamma, r);
        }
        out << (sum / N) << "\n";
    }
    out.close();
}

// -----------------------------------------------------------------------------
// main: punto di ingresso del programma
//   - Inizializza RNG
//   - Definisce parametri di simulazione
//   - Per ogni N in {1,2,10,100}:
//       genera file per Uniforme, Esponenziale e Cauchy
//   - Salva stato RNG e termina
// -----------------------------------------------------------------------------
int main() {
    // ----------- Inizializzazione del generatore di numeri casuali -----------
    Random rnd;
    int p1, p2;
    ReadPrimes("Primes", p1, p2);
    InitializeRandom(rnd, "seed.in", p1, p2);

    // -------------------- Parametri di simulazione -----------------------------
    const int M = 10000;            // numero di realizzazioni per ciascun N
    const int Ns[4] = {1, 2, 10, 100}; // valori di N da testare
    const double lambda = 1.0;      // parametro per distribuzione esponenziale
    const double mu = 0.0;          // parametro di posizione per Cauchy
    const double gamma = 1.0;       // parametro di scala per Cauchy

    // --------- Ciclo sui diversi valori di N per generare i dati -----------
    for (int idx = 0; idx < 4; ++idx) {
        int N = Ns[idx];
        // Costruzione dinamica dei nomi dei file
        ostringstream fStd, fExp, fCau;
        fStd << "Es01_2Standard_" << N << ".txt";
        fExp << "Es01_2Exp_"      << N << ".txt";
        fCau << "Es01_2Cauchy_"   << N << ".txt";

        // Generazione delle medie campionarie
        ComputeUniformAverages(rnd, M, N, fStd.str());
        ComputeExpAverages    (rnd, M, N, lambda, fExp.str());
        ComputeCauchyAverages (rnd, M, N, mu, gamma, fCau.str());

        // Conferma a video dell'avvenuta scrittura
        cout << "Dati generati per N=" << N
             << " nei file: " << fStd.str() << ","
             << fExp.str() << "," << fCau.str() << endl;
    }

    // Salvataggio dello stato del RNG per esecuzioni future
    rnd.SaveSeed();
    return 0;
}
