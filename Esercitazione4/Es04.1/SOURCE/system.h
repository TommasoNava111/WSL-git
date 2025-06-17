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
 * system.h: Dichiarazione della classe System                  *
 *                                                              *
 * Questo file contiene:                                        *
 *  - Inclusioni delle librerie necessarie (armadillo per        *
 *    operazioni vettoriali, iostream, fstream, ecc.).          *
 *  - Dichiarazione delle classi Particle e Random.             *
 *  - Definizione della classe System, che gestisce l’intero    *
 *    processo di simulazione (MD o MC).                        *
 *                                                              *
 * Struttura del file:                                           *
 *                                                              *
 * 1) Guardie di inclusione                                      *
 *    - #ifndef __System__, #define __System__, #endif           *
 *      per evitare inclusioni multiple.                        *
 *                                                              *
 * 2) Inclusioni                                                *
 *    - Librerie standard: iostream, iomanip, fstream, string,   *
 *      stdlib (per exit).                                       *
 *    - <armadillo>: per vettori e matrici (arma::vec, arma::field).*
 *    - "particle.h": dichiarazione della classe Particle, che    *
 *      rappresenta singole particelle (posizione, velocità,     *
 *      spin).                                                   *
 *    - "random.h": dichiarazione della classe Random, che fornisce*
 *      metodi per numeri casuali uniformi e gaussiani.           *
 *                                                              *
 * 3) Dichiarazione della classe System                          *
 *                                                              *
 * Attributi privati (campi dati):                               *
 *  - Costanti e flag generali:                                  *
 *      • const int _ndim = 3: numero di dimensioni (3D).        *
 *      • bool _restart: flag per restart da file precedente.     *
 *      • int _sim_type: tipo di simulazione (0=MD NVE, 1=MC NVT,  *
 *         2=Ising MRT², 3=Ising Gibbs).                           *
 *      • int _npart, _nblocks, _nsteps: numero di particelle,    *
 *         blocchi, step per blocco.                              *
 *      • int _nattempts, _naccepted: contatori di tentativi e    *
 *         accettazioni (per MC).                                  *
 *      • double _temp, _beta: temperatura e β=1/T.                *
 *      • double _rho, _volume: densità e volume.                  *
 *      • double _r_cut: cutoff per potenziale LJ.                 *
 *      • double _delta: passo di spostamento per MC.              *
 *      • double _J, _H: parametri del modello di Ising 1D.        *
 *  - Geometria                                                   *
 *      • arma::vec _side (3): lati della scatola (cubi).         *
 *      • arma::vec _halfside (3): metà dei lati (per pbc).         *
 *  - RNG                                                         *
 *      • Random _rnd: generatore di numeri casuali.               *
 *  - Stato delle particelle                                       *
 *      • arma::field<Particle> _particle: array di Particle,      *
 *         dimensione _npart.                                      *
 *      • arma::vec _fx, _fy, _fz: vettori di forze per ogni       *
 *         particella lungo x,y,z.                                 *
 *  - Proprietà da misurare                                        *
 *      • int _nprop: numero totale di osservabili (inclusi bins).*
 *      • bool flag per ogni proprietà: _measure_penergy,          *
 *        _measure_kenergy, _measure_tenergy (energie),            *
 *        _measure_temp (temperatura), _measure_pressure (P),      *
 *        _measure_gofr (radial distribution),                     *
 *        _measure_pofv (distribuzione velocità),                  *
 *        _measure_magnet (magnetizzazione),                        *
 *        _measure_cv (calore specifico), _measure_chi            *
 *        (suscettività).                                           *
 *      • int _index_*: indici di partenza nel vettore di misura    *
 *        (_measurement) per ogni proprietà.                         *
 *      • int _n_bins, _n_bins_v: numero di bin per gofr e pofv.     *
 *      • double _bin_size, _bin_size_v: ampiezza di ogni bin.       *
 *      • double _vtail, _ptail: correzioni di coda LJ (energia e   *
 *        pressione).                                                *
 *      • arma::vec _block_av, _global_av, _global_av2,              *
 *        _average, _measurement: vettori per blocchi e statistiche. *
 *      • int particelle: contatore temporaneo per hist pofv.        *
 *                                                              *
 * Metodi pubblici:                                              *
 *  - int get_nbl(): ritorna _nblocks.                             *
 *  - int get_nsteps(): ritorna _nsteps.                           *
 *  - void initialize(): legge input.dat, configura RNG,            *
 *    carica o genera posizioni e, se MD, velocità iniziali.        *
 *  - void initialize_properties(): legge properties.dat, crea      *
 *    file di output per ogni osservabile, alloca vettori,           *
 *    inizializza indici e flag.                                     *
 *  - void finalize(): salva configuration finale, seed RNG,         *
 *    scrive messaggi di chiusura su output.dat.                       *
 *  - void write_configuration(): salva configurazione corrente       *
 *    su ../OUTPUT/CONFIG/config.xyz (posizioni nuove e vecchie),       *
 *    o config.spin per Ising.                                            *
 *  - void write_XYZ(int nconf): salva snapshot istantaneo in file      *
 *    config_nconf.xyz (XYZ).                                               *
 *  - void read_configuration(): legge posizioni da ../INPUT/CONFIG/        *
 *    config.xyz, imposta _particle(i).pos nuove; se restart e Ising,      *
 *    legge spin da config.spin.                                               *
 *  - void initialize_velocities(): se restart, imposta velocità e      *
 *    pos “vecchie” da conf-1.xyz; altrimenti, genera velocità gaussiane,*
 *    rimuove drift, scalfa per T, imposta x_old = x(t) – v·Δt per Verlet. *
 *  - void step(): se sim_type==0 chiama Verlet(); altrimenti, per MC,      *
 *    sceglie particella casuale e chiama move(i).                               *
 *  - void block_reset(int blk): se blk>0 scrive “Block completed: blk” su      *
 *    output.dat, quindi azzera _block_av per nuovo blocco.                          *
 *  - void measure(): azzera _measurement;                                *
 *    • Se _measure_penergy/PRESSURE/goFR abilitati, calcola pot. LJ e         *
 *      viriale sommando su coppie i<j.                                             *
 *    • Se _measure_pofv abilitato, costruisce istogramma grezzo delle             *
 *      velocità in N_BINS_V bins su [0, 4·T], incrementando _measurement.          *
 *    • Se _measure_kenergy, calcola K = ½ Σ v_i².                                  *
 *    • Se _measure_tenergy, U+K o energia Ising.                                     *
 *    • Se _measure_temp, T = (2/3) K.                                                 *
 *    • Se _measure_pressure, P = ρ (2/3) K + [viriale + p_tail] / V.                 *
 *    • Somma _measurement in _block_av.                                              *
 *  - void averages(int blk):                                                                *
 *    • Calcola _average = _block_av / N_STEPS.                                              *
 *    • Aggiorna _global_av += _average, _global_av2 += _average².                           *
 *    • Per ogni proprietà abilitata, apre il relativo file (pot, kin, total, temp,         *
 *      pressure) e scrive:                                                                  *
 *        – blk (numero blocco), _average[index], media progressiva (_global_av[index]/blk), *
 *        – errore calcolato con error(_global_av[index], _global_av2[index], blk).            *
 *    • Se _measure_pofv e blk == NBLOCKS, apre pofv.dat e per ogni bin i:                 *
 *        – bin_center = (i+0.5)*_bin_size_v                                                *
 *        – avg = _global_av[_index_pofv+i]/blk                                             *
 *        – err = error(_global_av[_index_pofv+i], _global_av2[_index_pofv+i], blk)           *
 *        – Scrive bin_center, avg, err su pofv.dat.                                          *
 *    • Azzera visitatori di blocco con block_reset().                                       *
 *    • Aggiunge frazione di accettazione su acceptance.dat.                                   *
 *  - double error(double acc, double acc2, int blk):                                            *
 *    • Se blk ≤ 1, restituisce 0; altrimenti, calcola                                  *
 *      sqrt(|acc2/blk – (acc/blk)²| / blk).                                                 *
 *  - void move(int part): se sim_type=1 (MC LJ), propone spostamento casuale di particle[part];*
 *    chiama metro(part) per decidere accettazione e aggiorna stato di particle o ripristina. *
 *    Se sim_type indica Ising MRT²/ Gibbs, chiama metro(part) per flip spin.                     *
 *  - bool metro(int part): calcola ΔE per singola particella/spin (LJ o Ising) e accetta con    *
 *    probabilità exp(–β ΔE).                                                                   *
 *  - double pbc(double position, int i):                                *
 *    • Restituisce posizione ridotta nella scatola:                            *
 *      position – _side(i) * rint(position / _side(i)).                         *
 *  - int pbc(int i): per Ising, Se i ≥ _npart, restituisce i–_npart; se i<0, restituisce i+_npart.*
 *  - void Verlet():                                      *
 *    • Calcola forze su ogni particella con Force(i,dim).              *
 *    • Per ogni i, aggiorna x_new = 2 x(t) – x_old + F_i Δt²; imposta posizioni,  *
 *      aggiorna velocità v = (x_new – x_old)/(2 Δt) con pbc, chiama                      *
 *      particle.acceptmove() e setta new positions.                                         *
 *  - double Force(int i, int dim):                                *
 *    • Scorre j≠i, calcola distanza minima periodic, se dr<r_cut,        *
 *      aggiunge componente di forza LJ: distance(dim)*(48/dr¹⁴ – 24/dr⁸).               *
 *  - double Boltzmann(int i, bool xnew):                             *
 *    • Somma contributi LJ (1/r¹² – 1/r⁶) per coppia (i,j) con j≠i se dr<r_cut.         *
 *      Restituisce 4·energy_i per Metropolis.                                             *
 *  - Attributo pubblico: int particelle;                     *
 *    • Variabile temporanea usata in measure() per contare particelle con v<4T*.      *
 *                                                              *
 * 4) Fine della definizione della classe System e                    *
 *    chiusura di #endif per __System__.                                 *
 ****************************************************************/

#ifndef __System__
#define __System__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <armadillo>
#include <stdlib.h> //exit
#include "particle.h"
#include "random.h"

using namespace std;
using namespace arma;

class System {

private:
  const int _ndim = 3;  // Dimensionality of the system
  bool _restart;        // Flag indicating if the simulation is restarted
  int _sim_type;        // Type of simulation (e.g., Lennard-Jones, Ising)
  int _npart;           // Number of particles
  int _nblocks;         // Number of blocks for block averaging
  int _nsteps;          // Number of simulation steps in each block
  int _nattempts;       // Number of attempted moves
  int _naccepted;       // Number of accepted moves
  double _temp, _beta;  // Temperature and inverse temperature
  double _rho, _volume; // Density and volume of the system
  double _r_cut;        // Cutoff radius for pair interactions
  double _delta;        // Displacement step for particle moves
  double _J, _H;        // Parameters for the Ising Hamiltonian
  vec _side;            // Box dimensions
  vec _halfside;        // Half of box dimensions
  Random _rnd;          // Random number generator
  field <Particle> _particle; // Field of particle objects representing the system
  vec _fx, _fy, _fz;    // Forces on particles along x, y, and z directions

  // Properties
  int _nprop; // Number of properties being measured
  bool _measure_penergy, _measure_kenergy, _measure_tenergy;// Flags for measuring different energies
  bool _measure_temp, _measure_pressure, _measure_gofr;     // Flags for measuring temperature, pressure, and radial dist. function
  bool _measure_magnet, _measure_cv, _measure_chi;          // Flags for measuring magnetization, heat capacity, and susceptibility
  bool _measure_pofv;                                       // Flag for measuring the velocity modulus distribution
  int _index_penergy, _index_kenergy, _index_tenergy;       // Indices for accessing energy-related properties in vec _measurement
  int _index_temp, _index_pressure, _index_gofr;            // Indices for accessing temperature, pressure, and radial dist. function
  int _index_magnet, _index_cv, _index_chi;                 // Indices for accessing magnetization, heat capacity, and susceptibility
  int _index_pofv;                                          // Index for accessing velocity modulus distribution
  int _n_bins;           // Number of bins for radial distribution function
  int _n_bins_v;         // Number of bins for velocity modulus distribution
  double _bin_size;      // Size of bins for radial distribution function
  double _bin_size_v;    // Size of bins for velocity modulus distribution
  double _vtail, _ptail; // Tail corrections for energy and pressure
  vec _block_av;         // Block averages of properties
  vec _global_av;        // Global averages of properties
  vec _global_av2;       // Squared global averages of properties
  vec _average;          // Average values of properties
  vec _measurement;      // Measured values of properties

public: // Function declarations
  int get_nbl();              // Get the number of blocks
  int get_nsteps();           // Get the number of steps in each block
  void initialize();          // Initialize system properties
  void initialize_properties();// Initialize properties for measurement
  void finalize();            // Finalize system and clean up
  void write_configuration(); // Write final system configuration to XYZ file
  void write_XYZ(int nconf);  // Write system configuration in XYZ format on the fly
  void read_configuration();  // Read system configuration from file
  void initialize_velocities();// Initialize particle velocities
  void step();                // Perform a simulation step
  void block_reset(int blk);  // Reset block averages
  void measure();             // Measure properties of the system
  void averages(int blk);     // Compute averages of properties
  double error(double acc, double acc2, int blk); // Compute error
  void move(int part);        // Move a particle
  bool metro(int part);       // Perform Metropolis acceptance-rejection step
  double pbc(double position, int i); // Apply periodic boundary conditions for coordinates
  int pbc(int i);             // Apply periodic boundary conditions for spins
  void Verlet();              // Perform Verlet integration step
  double Force(int i, int dim); // Calculate force on a particle along a dimension
  double Boltzmann(int i, bool xnew); // Calculate Boltzmann factor for Metropolis acceptance
  int particelle;

};

#endif // __System__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
