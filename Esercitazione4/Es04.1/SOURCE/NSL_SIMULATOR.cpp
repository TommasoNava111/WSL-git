/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
/****************************************************************
 * main.cpp: Funzione principale per avviare la simulazione     *
 *                                                              *
 * In questo file si definisce il punto di ingresso dell'       *
 * applicazione. Viene creato un oggetto System, che gestisce   *
 * l'intera simulazione (MD o MC) per un gas di Lennard-Jones    *
 * o un modello di Ising 1D.                                      *
 *                                                              *
 * Struttura del file:                                           *
 *                                                              *
 * 1) Inclusioni                                                *
 *    - <iostream>: per operazioni di I/O standard.             *
 *    - "system.h": dichiarazione della classe System, che       *
 *      contiene tutti i metodi per l’inizializzazione,         *
 *      l’evoluzione temporale, le misure e le statistiche.      *
 *                                                              *
 * 2) Funzione main(int argc, char *argv[])                     *
 *    - Crea un contatore di configurazioni nconf=1 (per nome-   *
 *      file opzionali XYZ).                                     *
 *    - Istanzia il System (SYS).                                *
 *    - Chiama SYS.initialize(): legge i file di input,          *
 *      configura RNG e, se necessario, carica configurazioni   *
 *      precedenti.                                              *
 *    - Chiama SYS.initialize_properties(): apre i file di      *
 *      output per ogni osservabile (energia, temperatura,       *
 *      pressione, distribuzioni, ecc.) e inizializza i vettori  *
 *      interni per il data blocking.                             *
 *    - Chiama SYS.block_reset(0): azzera i contatori di blocco   *
 *      e segna l’inizio della simulazione.                       *
 *                                                              *
 * 3) Doppio ciclo di simulazione:                               *
 *    - Outer loop: for(i = 0; i < SYS.get_nbl(); i++)            *
 *        • Scorre i blocchi di simulazione (N_BLOCKS definito    *
 *          in properties o input.dat).                          *
 *    - Inner loop: for(j = 0; j < SYS.get_nsteps(); j++)         *
 *        • Esegue un singolo passo di simulazione:               *
 *            – SYS.step(): per MD chiama Verlet(); per MC        *
 *              effettua spostamenti/spin flips con Metropolis.  *
 *            – SYS.measure(): misura le proprietà (energia,      *
 *              pressione, istogrammi, ecc.) e accumula nel       *
 *              vettore _block_av.                                 *
 *            – Se j % 50 == 0: (commentato) SYS.write_XYZ(nconf)  *
 *              salvare la configurazione corrente in formato XYZ.*
 *              Incrementa nconf.                                  *
 *    - Dopo l’end del blocco (dopo N_STEPS):                      *
 *        • SYS.averages(i+1): calcola medie di blocco e valori    *
 *          progressivi per ogni proprietà e li scrive sui file.   *
 *        • SYS.block_reset(i+1): azzera _block_av e annota su     *
 *          output.dat la fine del blocco.                         *
 *                                                              *
 * 4) Al termine dei blocchi:                                    *
 *    - SYS.finalize(): salva la configurazione finale, il seed    *
 *      del RNG, e scrive le informazioni di chiusura su output.   *
 *    - return 0: termina il programma con successo.               *
 ****************************************************************/
#include <iostream>
#include "system.h"

using namespace std;

int main (int argc, char *argv[]){

  int nconf = 1;
  System SYS;
  SYS.initialize();
  SYS.initialize_properties();
  SYS.block_reset(0);

  for(int i=0; i < SYS.get_nbl(); i++){ //loop over blocks
    for(int j=0; j < SYS.get_nsteps(); j++){ //loop over steps in a block
      SYS.step();
      SYS.measure();
      if(j%50 == 0){
//        SYS.write_XYZ(nconf); //Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf++;
      }
    }
    SYS.averages(i+1);
    SYS.block_reset(i+1);
  }
  SYS.finalize();

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
