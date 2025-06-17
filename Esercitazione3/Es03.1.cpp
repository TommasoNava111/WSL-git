// Exercise 03.1: Plain Vanilla Option Pricing
// ---------------------------------------------------------------
// Teoria di Black-Scholes assume che il prezzo di un sottostante S(t)
// segua un moto browniano geometrico (GBM) con tasso privo di rischio r
// e volatilità σ costanti.
// La soluzione analitica per opzioni europee (Call/Put) a tempo t<T è:
//   C[S,t] = S N(d1) - K e^{-r(T-t)} N(d2)
//   P[S,t] = S (N(d1)-1) - K e^{-r(T-t)} (N(d2)-1)
// con d1 = [ln(S/K) + (r+σ^2/2)(T-t)]/(σ√(T-t)), d2 = d1 - σ√(T-t)
// e N(x) = (1/2)[1+erf(x/√2)] cumulativa di una N(0,1).
//
// Obiettivo Monte Carlo al tempo t=0 con parametri:
//   S0=100, K=100, T=1, r=0.1, σ=0.25
// due metodi di campionamento:
// 1) Sampling diretto di S(T) usando S(T)=S0 exp[(r-σ^2/2)T + σW(T)]
// 2) Sampling discretizzato sul percorso GBM: suddividere [0,T] in nSteps passi di dt = T/nSteps
//    e applicare ricorsione S(t+dt)=S(t) exp[(r-σ^2/2)dt + σ√dt Z]
//
// Per ciascun metodo e per call e put:
// - Simulare M estrazioni di S(T), payoff scontati exp(-rT)·max(±(S_T-K),0)
// - Organizzare in N blocks di L=M/N simulazioni ciascuno
// - In ogni blocco calcolare la media dei payoff (avgCall, avgPut)
// - Accumulare somme e somme quadratiche per calcolare la media cumulativa
//   e la deviazione standard della media
// - Salvare su quattro file:
//    direct_call.txt, direct_put.txt,
//    discretized_call.txt, discretized_put.txt
//   righe: [nBlock] [m - valore analitico] [errore]
//
// Questo approccio dimostra come la stima MC converge all'analitica
// e compara efficacia del campionamento diretto vs percorso discretizzato.
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;

// -----------------------------------------------------------------------------
// ReadPrimes: legge due numeri primi dal file indicato e li restituisce via
//             riferimento. Se l'apertura fallisce, stampa errore ed esce.
// -----------------------------------------------------------------------------
void ReadPrimes(const string& filename, int& p1, int& p2) {
    ifstream filePrimes(filename);
    if (!filePrimes.is_open()) {
        cerr << "ERRORE: impossibile aprire " << filename << endl;
        exit(1);
    }
    filePrimes >> p1 >> p2;
    filePrimes.close();
}

// -----------------------------------------------------------------------------
// InitializeRandom: legge il seme da seedFile, trova la keyword RANDOMSEED e
//                   quattro valori interi, quindi inizializza rnd con p1,p2.
// -----------------------------------------------------------------------------
void InitializeRandom(Random& rnd, const string& seedFile, int p1, int p2) {
    ifstream fileSeed(seedFile);
    if (!fileSeed.is_open()) {
        cerr << "ERRORE: impossibile aprire " << seedFile << endl;
        exit(1);
    }
    string key;
    int seedArr[4];
    while (fileSeed >> key) {
        if (key == "RANDOMSEED") {
            fileSeed >> seedArr[0] >> seedArr[1] >> seedArr[2] >> seedArr[3];
            rnd.SetRandom(seedArr, p1, p2);
            break;
        }
    }
    fileSeed.close();
}

// -----------------------------------------------------------------------------
// payoff_call: payoff di una call europea
// payoff_put : payoff di una put europea
// -----------------------------------------------------------------------------
inline double payoff_call(double S, double K) { return max(S - K, 0.0); }
inline double payoff_put (double S, double K) { return max(K - S, 0.0); }

// -----------------------------------------------------------------------------
// DirectOptionBlock:
//   - genera L simulazioni di prezzo terminale S_T con formula analitica
//   - calcola payoff scontati per call e put
//   - restituisce via riferimento media payoffs diretti
// -----------------------------------------------------------------------------
void DirectOptionBlock(Random& rnd, int L,
                       double S0, double K, double r, double sigma, double T,
                       double& avgCall, double& avgPut) {
    double sumCall = 0.0, sumPut = 0.0;
    for (int i = 0; i < L; ++i) {
        double Z = rnd.Gauss(0.0, 1.0);
        // Simulazione del prezzo terminale S_T
        double ST = S0 * exp((r - 0.5 * sigma * sigma) * T + sigma * sqrt(T) * Z);
        // Payoff scontato
        sumCall += exp(-r * T) * payoff_call(ST, K);
        sumPut  += exp(-r * T) * payoff_put (ST, K);
    }
    avgCall = sumCall / L;
    avgPut  = sumPut  / L;
}

// -----------------------------------------------------------------------------
// DiscreteOptionBlock:
//   - genera L cammini discretizzati in nSteps di dt = T/nSteps
//   - calcola payoff scontati per call e put basati sul prezzo finale
// -----------------------------------------------------------------------------
void DiscreteOptionBlock(Random& rnd, int L, int nSteps,
                          double S0, double K, double r, double sigma, double T,
                          double& avgCall, double& avgPut) {
    double sumCall = 0.0, sumPut = 0.0;
    double dt = T / double(nSteps);
    for (int i = 0; i < L; ++i) {
        double S = S0;
        for (int step = 0; step < nSteps; ++step) {
            double Z = rnd.Gauss(0.0, 1.0);
            S *= exp((r - 0.5 * sigma * sigma) * dt + sigma * sqrt(dt) * Z);
        }
        sumCall += exp(-r * T) * payoff_call(S, K);
        sumPut  += exp(-r * T) * payoff_put (S, K);
    }
    avgCall = sumCall / L;
    avgPut  = sumPut  / L;
}

// -----------------------------------------------------------------------------
// ComputeOptionEstimates:
//   - esegue N blocchi di simulazioni dirette e discretizzate
//   - calcola media cumulativa e incertezza per call e put
//   - scrive quattro file di output con le stime e l'errore
// -----------------------------------------------------------------------------
void ComputeOptionEstimates(Random& rnd,
                            int M, int N, int nSteps,
                            double S0, double K, double r, double sigma, double T) {
    int L = M / N;  // simulazioni per blocco
    // File di output
    ofstream fDirCall("direct_call.txt");
    ofstream fDirPut ("direct_put.txt");
    ofstream fDiscCall("discretized_call.txt");
    ofstream fDiscPut ("discretized_put.txt");
    if (!fDirCall || !fDirPut || !fDiscCall || !fDiscPut) {
        cerr << "ERRORE: impossibile aprire file di output" << endl;
        exit(1);
    }

    // Variabili per stima cumulativa
    double sumDirCall = 0.0, sumDirCall2 = 0.0;
    double sumDirPut  = 0.0, sumDirPut2  = 0.0;
    double sumDiscCall= 0.0, sumDiscCall2= 0.0;
    double sumDiscPut = 0.0, sumDiscPut2 = 0.0;

    for (int i = 0; i < N; ++i) {
        double avgDirCall, avgDirPut;
        double avgDiscCall, avgDiscPut;
        DirectOptionBlock(rnd, L, S0, K, r, sigma, T, avgDirCall, avgDirPut);
        DiscreteOptionBlock(rnd, L, nSteps, S0, K, r, sigma, T, avgDiscCall, avgDiscPut);

        // Aggiorno somme per stima cumulativa
        sumDirCall += avgDirCall;    sumDirCall2 += avgDirCall * avgDirCall;
        sumDirPut  += avgDirPut;     sumDirPut2  += avgDirPut  * avgDirPut;
        sumDiscCall+= avgDiscCall;   sumDiscCall2+= avgDiscCall * avgDiscCall;
        sumDiscPut += avgDiscPut;    sumDiscPut2 += avgDiscPut  * avgDiscPut;

        int n = i + 1;
        // medie e incertezze
        double mDirCall = sumDirCall / n;
        double eDirCall = (n>1) ? sqrt((sumDirCall2/n - mDirCall*mDirCall)/(n-1)) : 0.0;
        double mDirPut  = sumDirPut  / n;
        double eDirPut  = (n>1) ? sqrt((sumDirPut2 /n - mDirPut *mDirPut )/(n-1)) : 0.0;
        double mDiscCall= sumDiscCall/n;
        double eDiscCall= (n>1) ? sqrt((sumDiscCall2/n - mDiscCall*mDiscCall)/(n-1)) : 0.0;
        double mDiscPut = sumDiscPut /n;
        double eDiscPut = (n>1) ? sqrt((sumDiscPut2 /n - mDiscPut *mDiscPut )/(n-1)) : 0.0;

        // scrivo su file: (blocchi cumulativi) (media - valore esatto) errore
        // valori esatti presi da Black-Scholes: call≈14.9758, put≈5.45953
        fDirCall << n << " " << (mDirCall - 14.975790778311286) << " " << eDirCall  << "\n";
        fDirPut  << n << " " << (mDirPut  -  5.4595325819072364) << " " << eDirPut   << "\n";
        fDiscCall<< n << " " << (mDiscCall- 14.975790778311286) << " " << eDiscCall << "\n";
        fDiscPut << n << " " << (mDiscPut -  5.4595325819072364) << " " << eDiscPut  << "\n";
    }

    // chiudo file
    fDirCall.close(); fDirPut.close();
    fDiscCall.close(); fDiscPut.close();
}

// -----------------------------------------------------------------------------
// main: inizializza RNG, definisce parametri e chiama ComputeOptionEstimates
// -----------------------------------------------------------------------------
int main() {
    // Inizializzazione generatore casuale
    Random rnd;
    int p1, p2;
    ReadPrimes("Primes", p1, p2);
    InitializeRandom(rnd, "seed.in", p1, p2);

    // Parametri di simulazione
    const int M       = 100000;   // simulazioni totali
    const int Nblocks = 100;     // blocchi
    const int nSteps  = 100;     // passi per discretizzazione
    const double S0   = 100.0;   // prezzo iniziale sottostante
    const double K    = 100.0;   // strike
    const double T    = 1.0;     // maturità
    const double r    = 0.10;    // tasso privo di rischio
    const double sigma= 0.25;    // volatilità

    // Esecuzione calcolo prezzi opzioni
    ComputeOptionEstimates(rnd, M, Nblocks, nSteps, S0, K, r, sigma, T);

    // Salvataggio stato del generatore
    rnd.SaveSeed();
    return 0;
}
