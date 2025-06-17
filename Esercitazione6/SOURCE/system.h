/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __System__
#define __System__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <armadillo>
#include <stdlib.h> //exit
#include "particle.h"
#include <vector>
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
  double _sum_cv;      // Somma dei valori di Cv per blocco
  double _sum_cv2;     // Somma dei quadrati dei valori di Cv per blocco
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
  void set_T(double T);  // Imposta temperatura _temp e _beta = 1/T.
  void set_H(double H);  // Imposta campo magnetico esterno (Ising).
  void reset_measurements(); // Azzera _global_av e _global_av2 per iniziare nuova simulazione/serie di blocchi.
  int get_index_tenergy() const { return _index_tenergy; } // Indice energia totale per spin.
  int get_index_cv() const { return _index_cv; }  // Indice per E^2 (Cv).
  int get_index_chi() const { return _index_chi; } // Indice per suscettività
  int get_index_magnet() const { return _index_magnet; } // Indice per magnetizzazione.
  double get_global_av(int i) const { return _global_av(i); } // Accesso valore cumulativo di media.
  double get_global_av2(int i) const { return _global_av2(i); }  // Accesso valore cumulativo di media^2.
  double get_global_avc(int index, int blk) const {
    return _global_av(index) * blk;  // Restituisce somma cumulativa (non media) di proprietà 'index' su 'blk' blocchi
}
  int get_npart() const { return _npart; }  // Restituisce numero di particelle/spin.
  void set_sim_type(int sim_type) { _sim_type = sim_type; }  // Imposta tipo di simulazione.
  vector<double> _block_E_av;     // Accumula ⟨E⟩ per ogni blocco
  vector<double> _block_E2_av;
  double get_average(int index) const { return _average(index); } // Valore medio di blocco di proprietà 'index'.
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
