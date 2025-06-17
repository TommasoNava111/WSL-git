// Questo programma verifica l'uniformità del generatore Rannyu() tramite:
//
// 1) Stima di ⟨r⟩ = ∫₀¹ r dr = 1/2
//    • Funzione utilizzata: ComputeBlockStatistics(Random& rnd, int totalThrows, int nBlocks,
//         ofstream& outMean, ofstream& outVar)
//      – totalThrows (M): numero totale di estrazioni casuali.
//      – nBlocks (N): numero di blocchi per il metodo delle medie a blocchi.
//    • Variabili chiave all’interno:
//      – L = M / N: numero di lanci per ciascun blocco.
//      – cumMean, cumMean2: somme cumulative delle medie dei blocchi e dei loro quadrati,
//        per il calcolo di ⟨r⟩ e ⟨r⟩².
//      – blockMean = sum/L: media del blocco corrente.
//      – avgMean = cumMean/n: media cumulativa su n blocchi.
//      – errMean = sqrt((avgMean2 - avgMean²)/(n-1)): deviazione standard della media.
//    • Output su file Es01.1r.txt: colonne [numero lanci totali] [⟨r⟩ - 0.5] [errore].
//
// 2) Stima di σ² = ∫₀¹ (r - 1/2)² dr = 1/12
//    • Sempre in ComputeBlockStatistics:
//      – sumVar: somma di (r-0.5)² nel blocco.
//      – cumVar, cumVar2: somme cumulative delle medie di (r-0.5)² e dei loro quadrati.
//      – blockVar = sumVar/L: stima di σ² nel blocco.
//      – avgVar = cumVar/n: stima cumulativa di ⟨(r-0.5)²⟩.
//      – errVar = sqrt((avgVar2 - avgVar²)/(n-1)).
//    • Output su file Es01.1devi.txt: colonne [numero lanci totali] [⟨(r-0.5)²⟩ - 1/12] [errore].
//
// 3) Test di uniformità tramite χ² di Pearson
//    • Funzione: ComputeChiSquare(Random& rnd, int blockSize, int nBlocks,
//         int mBins, ofstream& outChi)
//      – blockSize: numero di estrazioni per blocco di test χ².
//      – nBlocks (chiBlocks): conteggio di blocchi su cui ripetere il test.
//      – mBins: numero di intervalli equispaziati su [0,1) per il conteggio.
//      – expected = blockSize / mBins: frequenza teorica in ciascun bin.
//    • Variabili interne:
//      – counts[mBins]: array intero per il conteggio delle occorrenze.
//      – chi2: somma di (counts[i] - expected)² / expected per i=0..mBins-1.
//    • Output su file Es01.1chi.txt: colonne [lanci finora] [valore χ²].
//
// Inizializzazione del generatore:
//   - ReadPrimes(string filename, int& p1, int& p2): legge due numeri primi per SetRandom().
//   - InitializeRandom(Random& rnd, string seedFile, int p1, int p2): trova "RANDOMSEED"
//     e quattro semi in seed.in, invocando rnd.SetRandom(seeds,p1,p2).
//
// Variabili globali in main():
//   M, N, chiBins, blockSize, chiBlocks: configurano l'estensione delle simulazioni.
//   outMean, outVar, outChi: stream di output per i tre test.

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;

// -----------------------------------------------------------------------------
// ReadPrimes: legge due numeri primi dal file indicato e li restituisce tramite
//             riferimento. Se il file non si apre, termina l'esecuzione col codice 1.
// -----------------------------------------------------------------------------
void ReadPrimes(const string& filename, int& p1, int& p2) {
    ifstream filePrimes(filename);
    if (!filePrimes.is_open()) {
        cerr << "ERRORE: Impossibile aprire il file " << filename << endl;
        exit(1);
    }
    // Il file contiene due interi separati da spazio o newline
    filePrimes >> p1 >> p2;
    filePrimes.close();
}

// -----------------------------------------------------------------------------
// InitializeRandom: legge il seme dal file di configurazione e inizializza
//                   il generatore Random con i primi p1 e p2 forniti.
// -----------------------------------------------------------------------------
void InitializeRandom(Random& rnd, const string& seedFile, int p1, int p2) {
    ifstream fileSeed(seedFile);
    if (!fileSeed.is_open()) {
        cerr << "ERRORE: Impossibile aprire il file " << seedFile << endl;
        exit(1);
    }
    string keyword;
    int seedArray[4];
    // Cerco la riga che inizia con "RANDOMSEED"
    while (fileSeed >> keyword) {
        if (keyword == "RANDOMSEED") {
            // Leggo quattro valori interi di seme
            fileSeed >> seedArray[0] >> seedArray[1] >> seedArray[2] >> seedArray[3];
            // Inizializzo il generatore
            rnd.SetRandom(seedArray, p1, p2);
            break;
        }
    }
    fileSeed.close();
}

// -----------------------------------------------------------------------------
// ComputeBlockStatistics:
//   - Genera M numeri casuali divisi in N blocchi
//   - Calcola la media di r e la media di (r-0.5)^2 in ogni blocco
//   - Effettua la media cumulativa e l'errore statistico (deviazione standard)
//   - Scrive su file i risultati per ogni blocco
// -----------------------------------------------------------------------------
void ComputeBlockStatistics(Random& rnd,
                            int totalThrows,    // numero totale di lanci M
                            int nBlocks,        // numero di blocchi N
                            ofstream& outMean,  // file per i risultati di <r>
                            ofstream& outVar) { // file per i risultati di <(r-0.5)^2>
    int L = totalThrows / nBlocks;  // lanci per blocco

    // Variabili per accumulare somme e somme al quadrato
    double cumMean = 0.0, cumMean2 = 0.0;
    double cumVar  = 0.0, cumVar2  = 0.0;

    // Ciclo sui blocchi
    for (int i = 0; i < nBlocks; ++i) {
        double sum      = 0.0;  // somma dei valori r
        double sumVar   = 0.0;  // somma dei quadrati delle deviazioni

        // Ciclo sui lanci all'interno del blocco
        for (int j = 0; j < L; ++j) {
            double r = rnd.Rannyu();      // numero casuale uniforme [0,1)
            sum += r;
            double dev = r - 0.5;        // deviazione da 0.5
            sumVar += dev * dev;         // quadrato della deviazione
        }

        // Media nel blocco di r e di (r-0.5)^2
        double blockMean = sum / L;
        double blockVar  = sumVar / L;

        // Aggiornamento somme cumulative
        cumMean += blockMean;
        cumMean2 += blockMean * blockMean;
        cumVar  += blockVar;
        cumVar2  += blockVar  * blockVar;

        int n = i + 1;  // numero di blocchi considerati finora

        // Media cumulativa di <r>
        double avgMean  = cumMean  / n;
        double avgMean2 = cumMean2 / n;
        // Errore statistico: se n>1 uso la deviazione standard campionaria
        double errMean = (n > 1)
            ? sqrt((avgMean2 - avgMean*avgMean) / (n - 1))
            : 0.0;

        // Media cumulativa di <(r-0.5)^2>
        double avgVar  = cumVar  / n;
        double avgVar2 = cumVar2 / n;
        double errVar = (n > 1)
            ? sqrt((avgVar2 - avgVar*avgVar) / (n - 1))
            : 0.0;

        // Scrittura su file:
        // outMean: (lanci totali finora) (〈r〉 - 0.5) errore
        outMean << n * L << " " << (avgMean - 0.5) << " " << errMean << "\n";
        // outVar: (lanci totali finora) (〈(r-0.5)^2〉 - 1/12) errore
        outVar  << n * L << " " << (avgVar - 1.0/12) << " " << errVar  << "\n";
    }
}

// -----------------------------------------------------------------------------
// ComputeChiSquare:
//   - Divide lanci in blocchi di dimensione blockSize
//   - Conta occorrenze in m bin equispaziati
//   - Calcola il test del chi quadro per ogni blocco
//   - Scrive su file il valore di chi2 cumulato
// -----------------------------------------------------------------------------
void ComputeChiSquare(Random& rnd,
                      int blockSize,    // lanci per blocco di chi2
                      int nBlocks,      // numero di blocchi
                      int mBins,        // numero di intervalli
                      ofstream& outChi) { // file per i risultati chi2
    double expected = double(blockSize) / mBins; // valore atteso in ciascun bin

    // Ciclo sui blocchi di test
    for (int i = 0; i < nBlocks; ++i) {
        int counts[100] = {0};  // array per conteggio (assume mBins<=100)

        // Genero blockSize numeri e conto il bin di appartenenza
        for (int j = 0; j < blockSize; ++j) {
            double r = rnd.Rannyu();
            int bin = int(r * mBins);
            counts[bin]++;
        }

        // Calcolo statistica chi2 per il blocco
        double chi2 = 0.0;
        for (int k = 0; k < mBins; ++k) {
            double diff = counts[k] - expected;
            chi2 += diff * diff / expected;
        }

        // Scrivo su file: (lanci totali finora) valore chi2
        int throwsSoFar = (i + 1) * blockSize;
        outChi << throwsSoFar << " " << chi2 << "\n";
    }
}

int main(int argc, char *argv[]) {
    // -------------------- Inizializzazione RNG --------------------
    Random rnd;         // istanza del generatore
    int p1, p2;         // primi per il SetRandom
    ReadPrimes("Primes", p1, p2);
    InitializeRandom(rnd, "seed.in", p1, p2);

    // -------------------- Parametri di simulazione --------------------
    const int M          = 100000; // numero totale di lanci per media e varianza
    const int N          = 100;    // numero di blocchi per media e varianza
    const int chiBins    = 100;    // numero di bin per chi quadro
    const int blockSize  = 10000;  // lanci per blocco per chi quadro
    const int chiBlocks  = 100;    // numero di blocchi per chi quadro

    // -------------------- Apertura file di output --------------------
    ofstream outMean("Es01.1r.txt");    // file per risultati di <r>
    ofstream outVar ("Es01.1devi.txt"); // file per risultati di varianza
    ofstream outChi ("Es01.1chi.txt");  // file per risultati chi2
    if (!outMean || !outVar || !outChi) {
        cerr << "ERRORE: Impossibile aprire uno o più file di output" << endl;
        return 1;
    }

    // -------------------- Esecuzione calcoli --------------------
    ComputeBlockStatistics(rnd, M, N, outMean, outVar);
    ComputeChiSquare      (rnd, blockSize, chiBlocks, chiBins, outChi);

    // -------------------- Salvataggio stato RNG e uscita --------------------
    rnd.SaveSeed();
    return 0;
}