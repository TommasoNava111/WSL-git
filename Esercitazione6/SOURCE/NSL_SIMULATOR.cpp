#include <iostream>
#include "system.h"
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

int main (int argc, char *argv[]){
  // Creo un array di temperature da 0.5 a 2.0 con passo 0.1
  vector<double> temperatures;
    for (double T = 0.5; T <= 2.0; T += 0.1) {
        temperatures.push_back(T);
    }
// Definisco due valori di campo esterno: H=0 (per U, C, X) e H=0.02 (per M)
    double H_values[] = {0.0, 0.02};
    // Definisco due stringhe per distinguere i casi "Gibbs" (_gibbs) e "Metropolis" (_metropolis)
    string methods[] = {"_gibbs", "_metropolis"};  // corrispondono a sim_type=3 e sim_type=2
// Ciclo prima sul metodo (Gibbs o Metropolis)
    for (string method : methods) {
      // Apro i file di output per U, C, X, M
        // Ogni file contiene: temperatura, valore_medio, errore_statistico
        ofstream U_file("../OUTPUT/U_vs_T" + method + ".dat");
        ofstream C_file("../OUTPUT/C_vs_T" + method + ".dat");
        ofstream X_file("../OUTPUT/X_vs_T" + method + ".dat");
        ofstream M_file("../OUTPUT/M_vs_T" + method + ".dat");
// Scrivo l'header nei file (commento iniziale)
        U_file << "# Temperature U U_error" << endl;
        C_file << "# Temperature C C_error" << endl;
        X_file << "# Temperature X X_error" << endl;
        M_file << "# Temperature M M_error" << endl;
// Ciclo su tutte le temperature
        for (double T : temperatures) {
            /// --- PARTE 1: simulazione con H = 0 per calcolare U, C, X ---

            // 1) Creo un oggetto System per il modello Ising
            System SYS;
            SYS.initialize(); // legge input.dat, inizializza RNG, carica configurazione iniziale
            SYS.set_T(T);  // imposto la temperatura corrente (aggiorna β = 1/T)
            SYS.set_H(0.0);  // imposto H=0 per energia, capacità e suscettività
             // Scelgo tra Gibbs e Metropolis in base alla stringa "method"
            if (method == "_gibbs") SYS.set_sim_type(3); // 3 = Gibbs
            else SYS.set_sim_type(2); // 2 = Metropolis
            SYS.initialize_properties(); // abilita le misure di energia, C, X, M, ecc.
            SYS.reset_measurements(); // azzera somme globali e di blocco
// Preparo vettori in cui salvare i risultati di ciascun blocco
            int nbl = SYS.get_nbl(); // numero di blocchi
            double beta = 1.0 / T;
            vector<double> C_blocks;  // capacità termica (media di blocco)
            vector<double> U_blocks; // energia interna per spin (media di blocco)
            vector<double> X_blocks;  // suscettività (media di blocco)
            vector<double> M_blocks; //magnetizzazione

 // Ciclo sui blocchi (Data Blocking)
for (int i = 0; i < SYS.get_nbl(); i++) {
  // In ciascun blocco, eseguo nsteps misure
    for (int j = 0; j < SYS.get_nsteps(); j++) {
        SYS.step(); // tenta un flip di spin (Gibbs o Metropolis)
        SYS.measure(); // misura H, H^2, M, M^2, ecc.
    }
    // Alla fine del blocco, aggiorno le medie di blocco e accumulate globali
    SYS.averages(i + 1);

    
   // Ora estraggo la media di energia interna 'per spin' calcolata dal blocco
    double U_block = SYS.get_average(SYS.get_index_tenergy());
    // SYS.get_index_tenergy() restituisce l'indice in _average di TOTAL_ENERGY per spin
    // Notare: in measure() TOTAL_ENERGY era già normalizzato per spin

    // E estraggo la media di H^2 (che in measure() avevo salvato in index_cv)
    double Etot2_block = SYS.get_average(SYS.get_index_cv());
    // Etot2_block = <(Energia totale)^2> (ma energia totale = per spin * N)

    // Calcolo la capacità termica di ciascun blocco secondo la formula:
    // C = β^2 [ <E^2> - <E>^2 ]   ma qui devo stare attento alle normalizzazioni.
    // Etot2_block = < (E_totale)^2 >. U_block = <E_per_spin>.
    // Dunque U_block * N = <E_totale>, e (U_block*N)^2 = <E_totale>^2.
    int N = SYS.get_npart();  // numero di spin
    double beta = 1.0 / T;
    double C_block = beta * beta * (Etot2_block - pow(U_block * N, 2)) / N;
     // Spiego i termini:
    // - Etot2_block è media di (E_totale)^2
    // - pow(U_block * N,2) è (media di E_totale)^2
    // - la differenza è Var(E_totale).
    // - Divido per N per scrivere C per spin.
    // Aggiungo i valori di blocco nei rispettivi vettori
    C_blocks.push_back(C_block);
    U_blocks.push_back(U_block);
    
     // Calcolo la suscettività di ciascun blocco:
    // SYS.get_index_chi() contiene < β ( ∑ s_i )^2 / N > calcolato in measure()
    double X_block = SYS.get_average(SYS.get_index_chi());
    X_blocks.push_back(X_block);
    // Resetto le somme di blocco per la prossima iterazione
    SYS.block_reset(i + 1);
}

 // Dopo aver raccolto i risultati di tutti i blocchi,
// calcolo media ed errore statistico per U, C, X:

// Energia U media su tutti i blocchi:
double U_avg = accumulate(U_blocks.begin(), U_blocks.end(), 0.0) / nbl;
// Deviazione std dei blocchi e divisione per sqrt(nbl)
double U_err = sqrt(inner_product(U_blocks.begin(), U_blocks.end(), U_blocks.begin(), 0.0) / nbl - pow(U_avg, 2)) / sqrt(nbl);
// Capacità termica C
double C_avg = accumulate(C_blocks.begin(), C_blocks.end(), 0.0) / nbl;
double C_err = sqrt(inner_product(C_blocks.begin(), C_blocks.end(), C_blocks.begin(), 0.0) / nbl - pow(C_avg, 2)) / sqrt(nbl);
// Suscettività X
double X_avg = accumulate(X_blocks.begin(), X_blocks.end(), 0.0) / nbl;
double X_err = sqrt(inner_product(X_blocks.begin(), X_blocks.end(), X_blocks.begin(), 0.0) / nbl - pow(X_avg, 2)) / sqrt(nbl);
 // Scrivo i risultati sul file corrispondente
U_file << T << " " << U_avg << " " << U_err << endl;
C_file << T << " " << C_avg << " " << C_err << endl;
X_file << T << " " << X_avg << " " << X_err << endl;

          
   // --- PARTE 2: simulazione con H = 0.02 per calcolare M ---

  // Ricreo un nuovo oggetto System (per evitare di mescolare con le misure precedenti)
  System SYS_M;
  SYS_M.initialize(); // inizializza RNG, legge input.dat, carica config iniziale
  SYS_M.set_T(T);  // imposto la stessa temperatura
  SYS_M.set_H(0.02); // imposto campo esterno h = 0.02
  if (method == "_gibbs") SYS_M.set_sim_type(3);
  else SYS_M.set_sim_type(2);
  SYS_M.initialize_properties(); // abilito le misure (tra cui magnetizzazione)
  // Equilibro la simulazione con 1000 step di "riscaldamento" (equilibrio)
  for(int i=0; i<1000; ++i) SYS_M.step();
   // Una volta equilibrato, posso resettare i contatori di blocco
  SYS_M.block_reset(0);
 // Preparo vettore in cui salvare la magnetizzazione media di ciascun blocco
 vector<double> M_values;
  // Ciclo su tutti i blocchi
  for (int i = 0; i < SYS_M.get_nbl(); i++) {
  for (int j = 0; j < SYS_M.get_nsteps(); j++) {
  SYS_M.step(); // muovo la configurazione
  SYS_M.measure(); // misuro M (magnetizzazione) per spin
    }
    SYS_M.averages(i + 1);
    // Estraggo la magnetizzazione media del blocco
    double M_block = SYS_M.get_average(SYS_M.get_index_magnet());
    M_values.push_back(M_block);
     // Resetto i contatori di blocco per la prossima iterazione
    SYS_M.block_reset(i + 1);
}
// Calcolo media ed errore tra i blocchi per M:
double M_avg = accumulate(M_values.begin(), M_values.end(), 0.0) / nbl;
double M_err = sqrt(inner_product(M_values.begin(), M_values.end(), M_values.begin(), 0.0) / nbl - pow(M_avg, 2)) / sqrt(nbl);
// Scrivo i risultati sul file di M
M_file << T << " " << M_avg << " " << M_err << endl;
        }

// Chiudo i file di output prima di passare oltre
        U_file.close();
        C_file.close();
        X_file.close();
        M_file.close();
    }
  
  return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
