// Es01.3.cpp
//
// Questo programma simula l'esperimento di Buffon per stimare il valore di π
// tramite il lancio casuale di un ago di lunghezza L su linee parallele distanti d.
// Il codice è organizzato come segue:
//
// 1) ReadPrimes(const string& filename, int& p1, int& p2)
//    - Legge due numeri primi dal file "Primes" per i parametri di SetRandom().
//    - Variabili in output: p1, p2 utilizzate in InitializeRandom().
//
// 2) InitializeRandom(Random& rnd, const string& seedFile, int p1, int p2)
//    - Cerca la keyword "RANDOMSEED" in seed.in e legge quattro semi interi.
//    - Inizializza rnd con rnd.SetRandom(seedArr,p1,p2) per assicurare la riproducibilità.
//
// 3) EstimatePiBlock(Random& rnd, double L, double d, int nThrows, ofstream& angoli)
//    - Parametri:
//      • L: lunghezza dell'ago.
//      • d: distanza tra le linee (d > L).
//      • nThrows: numero di tiri nel blocco.
//      • angoli: stream per salvare gli angoli θ generati.
//    - Per ogni lancio:
//      a) Genera x ∈ [0, d/2] con x = (d/2)*rnd.Rannyu().
//      b) Genera un angolo θ uniforme su [0, π/2] usando il metodo del cerchio
//         inscribendo (u,v) ∈ quadrato e calcolando sinθ = |u|/√(u²+v²).
//      c) Se x ≤ (L/2)*sinθ allora l'ago interseca la linea (colpo++).
//      d) Salva θ su angoli.txt tramite angoli << asin(sinθ) << endl;
//    - Output: stima di π nel blocco = (2*L*nThrows)/(colpi*d).
//
// 4) ComputePiEstimates(Random& rnd, double L, double d, int M, int N,
//                        const string& outFile, ofstream& angoli)
//    - Parametri:
//      • M: numero totale di tiri.
//      • N: numero di blocchi.
//      • outFile: nome del file Es01_3.txt per le stime cumulative.
//    - Divide M in N blocchi di nThrows = M/N tiri ciascuno.
//    - Per ciascun blocco invoca EstimatePiBlock e accumula sum e sum2:
//      • cumSum += π_block; cumSum2 += π_block².
//    - Calcola media cumulativa media = cumSum / nBlocksSoFar;
//      errore = √((media2 - media²)/(nBlocksSoFar-1)).
//    - Scrive su outFile: [tiri cumulativi] [media π] [errore].
//
// Variabili principali in main():
//   L = 1.0, d = 1.5, M = 1e6, N = 100.
//   angoli.txt: conserva angoli di ogni lancio.
//   Es01_3.txt: risultati di π vs numero di tiri.
//
// Tutti i file di output sono pronti per generare il plot di π e della sua
// incertezza in funzione del numero di lanci totali secondo il metodo dei blocchi.
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;

// -----------------------------------------------------------------------------
// ReadPrimes: legge due numeri primi dal file indicato e li restituisce via
//             riferimento. Se il file non si apre, stampa errore ed esce.
// -----------------------------------------------------------------------------
void ReadPrimes(const string& filename, int& p1, int& p2) {
    ifstream filePrimes(filename);
    if (!filePrimes.is_open()) {
        cerr << "ERRORE: Impossibile aprire il file " << filename << endl;
        exit(1);
    }
    filePrimes >> p1 >> p2;
    filePrimes.close();
}

// -----------------------------------------------------------------------------
// InitializeRandom: legge il seme dal file seedFile e inizializza il generatore
//                   Random con i primi p1 e p2. Cerca la keyword "RANDOMSEED"
//                   e poi quattro valori interi di seme.
// -----------------------------------------------------------------------------
void InitializeRandom(Random& rnd, const string& seedFile, int p1, int p2) {
    ifstream fileSeed(seedFile);
    if (!fileSeed.is_open()) {
        cerr << "ERRORE: Impossibile aprire il file " << seedFile << endl;
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
// EstimatePiBlock: esegue L lanci Monte Carlo per stimare pi in un singolo blocco
// Parametri:
//   rnd        - generatore di numeri casuali inizializzato
//   Lunghezza  - lunghezza dell'ago
//   d          - distanza tra linee parallele
//   L          - numero di lanci nel blocco
// Ritorna:
//   stima di pi basata sul numero di "colpi" (ago che tocca la linea)
// -----------------------------------------------------------------------------
double EstimatePiBlock(Random& rnd, double Lunghezza, double d, int L, ofstream& out) {
    int colpi = 0;
   // ofstream out("angoli");
    for (int i = 0; i < L; ++i) {
        // distanza dal centro dell'ago alla linea più vicina,
        // uniformemente distribuita tra 0 e d/2
        double x = (d / 2.0) * rnd.Rannyu();
        // generazione di un angolo uniforme sfruttando il cerchio unitario
        double u, v, r;
        do {
            u = rnd.Rannyu(-1.0, 1.0);
            v = rnd.Rannyu(-1.0, 1.0);
            r = u*u + v*v;
        } while (r > 1.0 || r == 0.0);
        // sin(theta) = |u|/sqrt(u^2+v^2)
        double sin_theta = fabs(u) / sqrt(r);
        out<< asin(sin_theta)<<endl;

        // proiezione orizzontale di metà ago
        double proiezione = (Lunghezza / 2.0) * sin_theta;
        // se la distanza x è inferiore o uguale alla proiezione, l'ago tocca la linea
        if (x <= proiezione) ++colpi;
        
    }
    // se non ci sono colpi, non ha senso calcolare pi
    if (colpi == 0) return 0.0;
    // formula di Buffon per stima di pi
    return (2.0 * Lunghezza * L) / (colpi * d);
}

// -----------------------------------------------------------------------------
// ComputePiEstimates: gestisce il calcolo a blocchi della stima di pi
//   - Divide M lanci in N blocchi di L = M/N lanci ciascuno
//   - Chiama EstimatePiBlock per ogni blocco
//   - Calcola media cumulativa e incertezza statistica a ogni passo
//   - Scrive i risultati su file di testo (lanci cumulativi, media, errore)
// -----------------------------------------------------------------------------
void ComputePiEstimates(Random& rnd,
                        double Lunghezza,
                        double d,
                        int M,
                        int N,
                        const string& filename, ofstream& angoli) {
    ofstream out(filename.c_str());
    if (!out.is_open()) {
        cerr << "ERRORE: Impossibile aprire il file " << filename << endl;
        exit(1);
    }

    int L = M / N; // lanci per blocco
    double cumSum = 0.0, cumSum2 = 0.0;

    for (int i = 0; i < N; ++i) {
        double piBlock = EstimatePiBlock(rnd, Lunghezza, d, L,angoli);
        cumSum  += piBlock;
        cumSum2 += piBlock * piBlock;

        int blocchiFinora = i + 1;
        double media   = cumSum  / blocchiFinora;
        double media2  = cumSum2 / blocchiFinora;
        double errore  = 0.0;
        if (blocchiFinora > 1) {
            errore = sqrt((media2 - media*media) / (blocchiFinora - 1));
        }
        out << (blocchiFinora * L) << " "  // numero di lanci cumulativi
            << media << " "                // stima media di pi
            << errore << "\n";           // incertezza statistica
    }
    out.close();
}

// -----------------------------------------------------------------------------
// main: punto di ingresso del programma
//   - Inizializza RNG tramite ReadPrimes e InitializeRandom
//   - Imposta parametri del problema
//   - Chiama ComputePiEstimates per generare file di output
//   - Salva stato del RNG e termina
// -----------------------------------------------------------------------------
int main() {
    // Inizializzazione del generatore di numeri casuali
    Random rnd;
    int p1, p2;
    ReadPrimes("Primes", p1, p2);
    InitializeRandom(rnd, "seed.in", p1, p2);

    // Parametri di simulazione
    const double Lunghezza = 1.0;   // lunghezza dell'ago
    const double d         = 1.5;   // distanza tra le linee
    const int M            = 1000000; // lanci totali
    const int N            = 100;    // blocchi
    ofstream angoli("angoli.txt");


    // File di output per la stima di pi
    const string filename = "Es01_3.txt";

    // Calcolo della stima di pi a blocchi e scrittura su file
    ComputePiEstimates(rnd, Lunghezza, d, M, N, filename,angoli);

    // Salvataggio dello stato del generatore per esecuzioni future
    rnd.SaveSeed();
    return 0;
}
