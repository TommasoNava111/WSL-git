// Es02.1.cpp
//
// Questo programma calcola via Monte Carlo l'integrale
//    I = ∫₀¹ (π/2) cos(π x / 2) dx = 1
// utilizzando due metodi di campionamento come richiesto:
//
// 1) Campionamento uniforme in [0,1]
//    • Funzione: ComputeBlockIntegral(Random& rnd, int L)
//      – Per ogni campione genera x = rnd.Rannyu() uniforme.
//      – Accumula sumUnif += f(x) con f(x) = (π/2) cos(π x / 2).
//    • Variabili di lavoro:
//      – L: numero di punti nel blocco (M/N).
//      – sumUnif, cumU, cumU2: somme e somme quadratiche per media e varianza.
//    • Uscita: nel file Es02.1unif.txt, per ogni blocco scrive
//      [campioni cumulativi] [⟨I⟩ - 1] [errore statistical]
//
// 2) Importance sampling con densità p(x) = 2(1 - x)
//    • Metodo di inversione: xIS = 1 - sqrt(1 - r), r = rnd.Rannyu().
//    • Peso campionamento: peso = f(xIS) / p(xIS).
//    • Accumula sumIS += peso.
//    • Variabili di lavoro:
//      – sumIS, cumI, cumI2: somme e somme quadratiche per media e varianza.
//    • Uscita: nel file Es02.1nonunif.txt, per ogni blocco scrive
//      [campioni cumulativi] [⟨I⟩ - 1] [errore statistical]
//
// 3) Funzione ComputeIntegralBlocks:
//    • Parametri: M (campioni totali), N (blocchi), rnd, nomi file.
//    • L = M/N: campioni per blocco.
//    • Per ogni blocco chiama ComputeBlockIntegral, aggiorna cumU, cumU2, cumI, cumI2.
//    • Calcola le medie cumulative meanU, meanI e gli errori errU, errI.
//    • Scrive su file i risultati per tracciare i plot richiesti.
//
// Variabili principali in main():
//    M = 1e4, N = 1e2, fileU = "Es02.1unif.txt", fileI = "Es02.1nonunif.txt".
//    Il generatore Random rnd si inizializza con ReadPrimes e InitializeRandom
//    per garantire la riproducibilità delle estrazioni.
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;

// -----------------------------------------------------------------------------
// ReadPrimes: legge due numeri primi dal file indicato e li restituisce via
//             riferimento. Se il file non si apre, stampa un errore ed esce.
// -----------------------------------------------------------------------------
void ReadPrimes(const string& filename, int& p1, int& p2) {
    ifstream filePrimes(filename);
    if (!filePrimes.is_open()) {
        cerr << "ERRORE: Impossibile aprire " << filename << endl;
        exit(1);
    }
    filePrimes >> p1 >> p2;
    filePrimes.close();
}

// -----------------------------------------------------------------------------
// InitializeRandom: legge il seme dal file seedFile, cerca la keyword
//                   "RANDOMSEED" e quattro interi, quindi inizializza rnd.
// -----------------------------------------------------------------------------
void InitializeRandom(Random& rnd, const string& seedFile, int p1, int p2) {
    ifstream fileSeed(seedFile);
    if (!fileSeed.is_open()) {
        cerr << "ERRORE: Impossibile aprire " << seedFile << endl;
        exit(1);
    }
    string keyword;
    int seedArr[4];
    while (fileSeed >> keyword) {
        if (keyword == "RANDOMSEED") {
            fileSeed >> seedArr[0] >> seedArr[1] >> seedArr[2] >> seedArr[3];
            rnd.SetRandom(seedArr, p1, p2);
            break;
        }
    }
    fileSeed.close();
}

// -----------------------------------------------------------------------------
// ComputeBlockIntegral: calcola l'integrale di f(x)=cos(pi x/2) su [0,1]
//   con campionamento uniforme e non-uniforme in un singolo blocco di L punti.
//   Ritorna la coppia {I_uniforme, I_importance}
//   - rnd      : generatore casuale
//   - L        : numero di campioni nel blocco
// -----------------------------------------------------------------------------
pair<double,double> ComputeBlockIntegral(Random& rnd, int L) {
    double sumUnif = 0.0;
    double sumIS   = 0.0;
    for (int i = 0; i < L; ++i) {
        // Campionamento uniforme su [0,1]
        double x = rnd.Rannyu();
        sumUnif += (M_PI/2.0) * cos(M_PI * x / 2.0);
        // Importance sampling con densità p(r) = 2*(1-r)
        double r = rnd.Rannyu();
        double xIS = 1.0 - sqrt(1.0 - r);
        double peso = (M_PI/2.0) * cos(M_PI * xIS / 2.0) / (2.0*(1.0 - xIS));
        sumIS += peso;
    }
    return { sumUnif / L, sumIS / L };
}

// -----------------------------------------------------------------------------
// ComputeIntegralBlocks: esegue N blocchi di M lanci per il calcolo dell'integrale
//   - rnd           : generatore casuale
//   - M             : numero totale di campioni
//   - N             : numero di blocchi
//   - filenameUnif  : nome file output per campionamento uniforme
//   - filenameIS    : nome file output per importance sampling
// -----------------------------------------------------------------------------
void ComputeIntegralBlocks(Random& rnd,
                            int M,
                            int N,
                            const string& filenameUnif,
                            const string& filenameIS) {
    ofstream outU(filenameUnif.c_str());
    ofstream outI(filenameIS.c_str());
    if (!outU || !outI) {
        cerr << "ERRORE: Impossibile aprire i file di output" << endl;
        exit(1);
    }

    int L = M / N; // campioni per blocco
    double cumU = 0.0, cumU2 = 0.0;
    double cumI = 0.0, cumI2 = 0.0;

    for (int i = 0; i < N; ++i) {
        auto [Iu, Ii] = ComputeBlockIntegral(rnd, L);
        cumU  += Iu;    cumU2 += Iu * Iu;
        cumI  += Ii;    cumI2 += Ii * Ii;

        int n = i + 1;
        double meanU = cumU / n;
        double meanU2 = cumU2 / n;
        double errU = (n>1)
            ? sqrt((meanU2 - meanU*meanU)/(n-1))
            : 0.0;

        double meanI = cumI / n;
        double meanI2 = cumI2 / n;
        double errI = (n>1)
            ? sqrt((meanI2 - meanI*meanI)/(n-1))
            : 0.0;

        // Scrittura: (campioni cumulativi) (mean - valore esatto=1) errore
        outU << n*L << " " << (meanU - 1.0) << " " << errU << "\n";
        outI << n*L << " " << (meanI - 1.0) << " " << errI << "\n";
    }
    outU.close();
    outI.close();
}

// -----------------------------------------------------------------------------
// main: inizializza RNG, definisce parametri e lancia il calcolo a blocchi
// -----------------------------------------------------------------------------
int main() {
    Random rnd;
    int p1,p2;
    ReadPrimes("Primes", p1, p2);
    InitializeRandom(rnd, "seed.in", p1, p2);

    const int M = 10000;    // campioni totali
    const int N = 100;      // blocchi
    const string fileU = "Es02.1unif.txt";
    const string fileI = "Es02.1nonunif.txt";

    ComputeIntegralBlocks(rnd, M, N, fileU, fileI);

    rnd.SaveSeed();
    return 0;
}
