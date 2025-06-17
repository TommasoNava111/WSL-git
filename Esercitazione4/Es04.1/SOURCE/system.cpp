/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
_/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
 
#include <cmath>
#include <cstdlib>
#include <string>
#include "system.h"
 
using namespace std;
using namespace arma;
 
/***************************************************************
 * System.cpp: Implementazione di simulazioni di Dynamics e     *
 * Monte Carlo per gas di Lennard-Jones e Modello di Ising 1D   *
 *                                                              *
 * In questo file sono definite le funzioni principali per:     *
 *                                                              *
 * 1) Inizializzazione del sistema                               *
 *    - System::initialize():                                   
 *      • Legge i parametri da input.dat (tipo di simulazione,  *
 *        temperatura, numero di particelle, densità, cutoff,   *
 *        delta, numero di blocchi, numero di step, ecc.)       *
 *      • Inizializza RNG e, se richiesto, carica una            *
 *        configurazione precedente da file; altrimenti genera   *
 *        posizioni iniziali e velocità (per MD).               *
 *    - System::initialize_velocities():                        
 *      • Se è un restart, legge velocità e posizioni “vecchie”  *
 *        da conf-1.xyz; altrimenti, assegna velocità casuali    *
 *        secondo una Gaussiana centrata su 0 e normalizza la    *
 *        somma delle velocità per avere energia cinetica risp.  *
 *        a T desiderata.                                         *
 *    - System::initialize_properties():                        
 *      • Legge properties.dat per capire quali variabili       
 *        estrarre (energia potenziale, cinetica, totale,         
 *        temperatura, pressione, radial distribution,           *
 *        distribuzione velocità, magnetizzazione,            *
 *        specific heat, suscettività).                         
 *      • Crea i file di output per ogni proprietà e alloca     *
 *        gli array interni (_measurement, _average, _block_av, *
 *        _global_av, _global_av2) con dimensione pari al numero *
 *        di osservabili (più eventuali bins per gofr, pofv, ecc.).*
 *                                                              *
 * 2) Evoluzione del sistema                                      *
 *    - System::step():                                          
 *      • Se sim_type == 0, chiama Verlet() per un passo MD;     
 *      • Altrimenti, per MC, sceglie a caso una particella e la *
 *        muove o le gira lo spin.                                *
 *      • Incrementa _nattempts (per contare proposte di MC).    *
 *    - System::Verlet():                                        
 *      • Calcola le forze per ogni particella tramite Force().  
 *      • Applica lo schema di integrazione di Verlet:           
 *        x_new = 2·x(t) – x(t–Δt) + F·Δt²;                      
 *      • Aggiorna velocità come (x_new – x_old)/(2Δt) e impone   *
 *        periodicità con pbc().                                   *
 *    - System::Force(int i, int dim):                           
 *      • Scorre tutte le coppie (i, j) con j ≠ i; calcola la     *
 *        distanza minima con condizioni periodiche;              
 *      • Se dr < r_cut, somma al gradiente LJ:                    
 *        F_dim += Δ_dim·(48/dr¹⁴ – 24/dr⁸).                       
 *    - System::move(int i):                                      
 *      • Proprone uno spostamento casuale (shift) in [-Δ, +Δ]^3  
 *        per MC; accetta o rifiuta tramite metro(i);             
 *      • Se sim_type indica Ising, chiama metro(i) per cambiare *
 *        lo spin.                                                 *
 *    - System::metro(int i):                                    
 *      • Calcola ΔE locale; accetta con probabilità exp(–βΔE).  
 *    - System::Boltzmann(int i, bool xnew):                      
 *      • Calcola l’energia LJ (pre-tagliata) associata alla       *
 *        particella i (data la posizione nuova o attuale).       *
 *    - System::pbc(double pos, int dim) / pbc(int i):           
 *      • Applica condizioni al contorno periodiche per posizioni *
 *        reali (riporta dentro la scatola) o indici Ising.       *
 *                                                              *
 * 3) Monitoraggio delle proprietà e Data Blocking                *
 *    - System::measure():                                       
 *      • Resetta _measurement a zero.                           
 *      • Se richiesto, calcola:                                 *
 *          – Energia potenziale sommando contributi LJ su tutte  *
 *            le coppie (i<j).                                    *
 *          – Velocità ≔ ||v_i|| per ogni particella; se abilitato,*
 *            costruisce un istogramma grezzo per p(v) su N_bins. *
 *          – Energia cinetica = ½ Σ m v_i².                      
 *          – Energia totale = K + U.                              *
 *          – Temperatura = (2/3) K (per MD).                     *
 *          – Pressione = ρ(2/3 K) + termine da viriale + tail.    *
 *      • Aggiorna i conteggi di istogrammi (goFR, pofv) con       *
 *        normalizzazione provvisoria all’interno di ciascun blocco.*
 *      • Somma il vettore _measurement in _block_av per blocco.   *
 *    - System::block_reset(int blk):                            
 *      • Se blk>0, stampa su output.dat “Block completed: blk”.  
 *      • Azzera _block_av per iniziare nuovo blocco.              
 *    - System::averages(int blk):                                
 *      • Calcola la media di blocco: _average = _block_av/N_steps. *
 *      • Aggiorna _global_av e _global_av2 con _average e _average².*
 *      • Per ogni proprietà abilitata, scrive su file:           
 *          – media di blocco (ACTUAL)                                 *
 *          – media progressiva (AVE)                                 *
 *          – errore statistico calcolato con la funzione error().   *
 *      • Per p(v) (POFV), agisce solo all’ultimo blocco: calcola  *
 *        centro bin, media progressiva e incertezza per ogni bin.  *
 *      • Aggiorna acceptance fraction in acceptance.dat.           *
 *    - System::error(double sum, double sum2, int blk):          
 *      • Se blk≤1 ritorna 0; altrimenti, σ = sqrt(|sum2/blk – (sum/blk)²|/blk).*
 *                                                              *
 * 4) Output di configurazioni                                       *
 *    - System::write_configuration():                            
 *      • Scrive su ../OUTPUT/CONFIG/config.xyz la configurazione  *
 *        corrente (posizioni “nuove”) e config/conf-1.xyz         *
 *        (posizioni “vecchie”) in unità ridotte. Se Ising, salva *
 *        spin sequence su config.spin.                            *
 *    - System::write_XYZ(int nconf):                              *
 *      • Salva la configurazione corrente in file config_nconf.xyz.*
 *    - System::read_configuration():                              *
 *      • Legge posizioni da ../INPUT/CONFIG/config.xyz e le imposta*
 *        come nuove; se restart e Ising, legge anche config.spin. *
 *                                                              *
 * 5) Utility                                                        *
 *    - System::finalize():                                        *
 *      • Scrive la configurazione finale e salva il seed RNG.    
 *    - System::get_nbl(): ritorna _nblocks.                      
 *    - System::get_nsteps(): ritorna _nsteps.                    
 ***************************************************************/

 
void System :: step(){ // Perform a simulation step
  if(_sim_type == 0) this->Verlet();  // Perform a MD step
  else for(int i=0; i<_npart; i++) this->move(int(_rnd.Rannyu()*_npart)); // Perform a MC step on a randomly choosen particle
  _nattempts += _npart; //update number of attempts performed on the system
  return;
}
 
void System :: Verlet(){
  double xnew, ynew, znew;
  for(int i=0; i<_npart; i++){ //Force acting on particle i
    _fx(i) = this->Force(i,0);
    _fy(i) = this->Force(i,1);
    _fz(i) = this->Force(i,2);
  }
  for(int i=0; i<_npart; i++){ //Verlet integration scheme (for the whole gas)
    xnew = this->pbc( 2.0 * _particle(i).getposition(0,true) - _particle(i).getposition(0,false) + _fx(i) * pow(_delta,2), 0);
    ynew = this->pbc( 2.0 * _particle(i).getposition(1,true) - _particle(i).getposition(1,false) + _fy(i) * pow(_delta,2), 1);
    znew = this->pbc( 2.0 * _particle(i).getposition(2,true) - _particle(i).getposition(2,false) + _fz(i) * pow(_delta,2), 2);
    _particle(i).setvelocity(0, this->pbc(xnew - _particle(i).getposition(0,false), 0)/(2.0 * _delta));
    _particle(i).setvelocity(1, this->pbc(ynew - _particle(i).getposition(1,false), 1)/(2.0 * _delta));
    _particle(i).setvelocity(2, this->pbc(znew - _particle(i).getposition(2,false), 2)/(2.0 * _delta));
    _particle(i).acceptmove(); // xold = xnew
    _particle(i).setposition(0, xnew);
    _particle(i).setposition(1, ynew);
    _particle(i).setposition(2, znew);
  }
  _naccepted += _npart;
  return;
}
 
double System :: Force(int i, int dim){
  double f=0.0, dr;
  vec distance;
  distance.resize(_ndim);
  for (int j=0; j<_npart; j++){
    if(i != j){
      distance(0) = this->pbc( _particle(i).getposition(0,true) - _particle(j).getposition(0,true), 0);
      distance(1) = this->pbc( _particle(i).getposition(1,true) - _particle(j).getposition(1,true), 1);
      distance(2) = this->pbc( _particle(i).getposition(2,true) - _particle(j).getposition(2,true), 2);
      dr = sqrt( dot(distance,distance) );
      if(dr < _r_cut){
        f += distance(dim) * (48.0/pow(dr,14) - 24.0/pow(dr,8)); //gradient formula for lennard-jones potential
      }
    }
  }
  return f;
}
 
void System :: move(int i){ // Propose a MC move for particle i
  if(_sim_type == 3){ //Gibbs sampler for Ising
    // TO BE FIXED IN EXERCISE 6
  } else {           // M(RT)^2
    if(_sim_type == 1){       // LJ system
      vec shift(_ndim);       // Will store the proposed translation
      for(int j=0; j<_ndim; j++){
        shift(j) = _rnd.Rannyu(-1.0,1.0) * _delta; // uniform distribution in [-_delta;_delta)
      }
      _particle(i).translate(shift, _side);  //Call the function Particle::translate
      if(this->metro(i)){ //Metropolis acceptance evaluation
        _particle(i).acceptmove();
        _naccepted++;
      } else _particle(i).moveback(); //If translation is rejected, restore the old configuration
    } else {                  // Ising 1D
      if(this->metro(i)){     //Metropolis acceptance evaluation for a spin flip involving spin i
        _particle(i).flip();  //If accepted, the spin i is flipped
        _naccepted++;
      }
    }
  }
  return;
}
 
bool System :: metro(int i){ // Metropolis algorithm
  bool decision = false;
  double delta_E, acceptance;
  if(_sim_type == 1) delta_E = this->Boltzmann(i,true) - this->Boltzmann(i,false);
  else delta_E = 2.0 * _particle(i).getspin() *
                 ( _J * (_particle(this->pbc(i-1)).getspin() + _particle(this->pbc(i+1)).getspin() ) + _H );
  acceptance = exp(-_beta*delta_E);
  if(_rnd.Rannyu() < acceptance ) decision = true; //Metropolis acceptance step
  return decision;
}
 
double System :: Boltzmann(int i, bool xnew){
  double energy_i=0.0;
  double dx, dy, dz, dr;
  for (int j=0; j<_npart; j++){
    if(j != i){
      dx = this->pbc(_particle(i).getposition(0,xnew) - _particle(j).getposition(0,1), 0);
      dy = this->pbc(_particle(i).getposition(1,xnew) - _particle(j).getposition(1,1), 1);
      dz = this->pbc(_particle(i).getposition(2,xnew) - _particle(j).getposition(2,1), 2);
      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);
      if(dr < _r_cut){
        energy_i += 1.0/pow(dr,12) - 1.0/pow(dr,6);
      }
    }
  }
  return 4.0 * energy_i;
}
 
double System :: pbc(double position, int i){ // Enforce periodic boundary conditions
  return position - _side(i) * rint(position / _side(i));
}
 
int System :: pbc(int i){ // Enforce periodic boundary conditions for spins
  if(i >= _npart) i = i - _npart;
  else if(i < 0)  i = i + _npart;
  return i;
}
 
void System :: initialize(){ // Initialize the System object according to the content of the input files in the ../INPUT/ directory
 
  int p1, p2; // Read from ../INPUT/Primes a pair of numbers to be used to initialize the RNG
  ifstream Primes("../INPUT/Primes");
  Primes >> p1 >> p2 ;
  Primes.close();
  int seed[4]; // Read the seed of the RNG
  ifstream Seed("../INPUT/seed.in");
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  _rnd.SetRandom(seed,p1,p2);
 
  ofstream couta("../OUTPUT/acceptance.dat"); // Set the heading line in file ../OUTPUT/acceptance.dat
  couta << "#   N_BLOCK:  ACCEPTANCE:" << endl;
  couta.close();
 
  ifstream input("../INPUT/input.dat"); // Start reading ../INPUT/input.dat
  ofstream coutf;
  coutf.open("../OUTPUT/output.dat");
  string property;
  double delta;
  while ( !input.eof() ){
    input >> property;
    if( property == "SIMULATION_TYPE" ){ //definisco il tipo di simulazione del sistema
      input >> _sim_type;
      if(_sim_type > 1){
        input >> _J;
        input >> _H;
      }
      if(_sim_type > 3){
        cerr << "PROBLEM: unknown simulation type" << endl;
        exit(EXIT_FAILURE);
      }
      if(_sim_type == 0)      coutf << "LJ MOLECULAR DYNAMICS (NVE) SIMULATION"  << endl;
      else if(_sim_type == 1) coutf << "LJ MONTE CARLO (NVT) SIMULATION"         << endl;
      else if(_sim_type == 2) coutf << "ISING 1D MONTE CARLO (MRT^2) SIMULATION" << endl;
      else if(_sim_type == 3) coutf << "ISING 1D MONTE CARLO (GIBBS) SIMULATION" << endl;
    } else if( property == "RESTART" ){ //il sistema è alla prima lettura o è gia stato letto una volta
      input >> _restart;
    } else if( property == "TEMP" ){
      input >> _temp;
      _beta = 1.0/_temp;
      coutf << "TEMPERATURE= " << _temp << endl;
    } else if( property == "NPART" ){
      input >> _npart;
      _fx.resize(_npart); //resetto il vettore forze lungo x con tutte le forze lungo x delle particelle
      _fy.resize(_npart); //resetto il vettore forze lungo x con tutte le forze lungo x delle particelle
      _fz.resize(_npart); //resetto il vettore forze lungo x con tutte le forze lungo x delle particelle
      _particle.set_size(_npart); //resetta il field delle particelle con il numero di particelle in input
      for(int i=0; i<_npart; i++){   //inizialiazzo tutte le particelle nel field
        _particle(i).initialize();
        if(_rnd.Rannyu() > 0.5) _particle(i).flip(); // to randomize the spin configuration
      }
      coutf << "NPART= " << _npart << endl;
    } else if( property == "RHO" ){
      input >> _rho;
      _volume = _npart/_rho; //il volume è la il numero di particella sulla densità
      _side.resize(_ndim); //setto il vettore dei lati come il numero di dimensioni
      _halfside.resize(_ndim);  //setto il vettore dei mezzi lati come il numero di dimensioin
      double side = pow(_volume, 1.0/3.0); //variabile lato del sistema cubico
      for(int i=0; i<_ndim; i++) _side(i) = side; //riempio il vettore dei lati
      _halfside=0.5*_side; //riempio il vettore dei mezzi lati
      coutf << "SIDE= ";
      for(int i=0; i<_ndim; i++){
        coutf << setw(12) << _side[i];
      }
      coutf << endl;
    } else if( property == "R_CUT" ){
      input >> _r_cut;
      coutf << "R_CUT= " << _r_cut << endl;
    } else if( property == "DELTA" ){
      input >> delta;
      coutf << "DELTA= " << delta << endl;
      _delta = delta;
    } else if( property == "NBLOCKS" ){
      input >> _nblocks;
      coutf << "NBLOCKS= " << _nblocks << endl;
    } else if( property == "NSTEPS" ){
      input >> _nsteps;
      coutf << "NSTEPS= " << _nsteps << endl;
    } else if( property == "ENDINPUT" ){
      coutf << "Reading input completed!" << endl;
      break;
    } else cerr << "PROBLEM: unknown input" << endl;
  }
  input.close();
  this->read_configuration();
  if(_sim_type==0) this->initialize_velocities();
  coutf << "System initialized!" << endl;
  coutf.close();
  return;
}
 
void System :: initialize_velocities(){
  double xold, yold, zold;
  if(_restart){
    ifstream cinf;
    cinf.open("../INPUT/CONFIG/conf-1.xyz");
    if(cinf.is_open()){
      string comment;
      string particle;
      int ncoord;
      cinf >> ncoord;
      if (ncoord != _npart){
        cerr << "PROBLEM: conflicting number of coordinates in input.dat & config.xyz not match!" << endl;
        exit(EXIT_FAILURE);
      }
      cinf >> comment;
      for(int i=0; i<_npart; i++){
        cinf >> particle >> xold >> yold >> zold; // units of coordinates in conf.xyz is _side
        _particle(i).setpositold(0, this->pbc(_side(0)*xold, 0));
        _particle(i).setpositold(1, this->pbc(_side(1)*yold, 1));
        _particle(i).setpositold(2, this->pbc(_side(2)*zold, 2));
      }
    } else cerr << "PROBLEM: Unable to open INPUT file conf-1.xyz"<< endl;
    cinf.close();
  } else {
    vec vx(_npart), vy(_npart), vz(_npart); //creo un vettore per dimensione lunghi pari al numero di particelle del sistema
    vec sumv(_ndim); //vettore lungo il numero di dimensioni del sistema -- ogni cella somma le velocità per dimensione delle particelle
    sumv.zeros(); //setto il vettore delle somme delle velocità a zero
    for (int i=0; i<_npart; i++){ //riempio ogni particella con le velocità nelle tre dimensioni
      vx(i) = _rnd.Gauss(0.,sqrt(_temp)); //generate velocities distributed as maxwell-boltzmann distribution
      vy(i) = _rnd.Gauss(0.,sqrt(_temp));
      vz(i) = _rnd.Gauss(0.,sqrt(_temp));
      sumv(0) += vx(i); //drift velocities of the gas -- averege velocity per particle
      sumv(1) += vy(i);
      sumv(2) += vz(i);
    }
    for (int idim=0; idim<_ndim; idim++) sumv(idim) = sumv(idim)/double(_npart); //velcoità per particella per ogni dimensione
    double sumv2 = 0.0, scalef;
    for (int i=0; i<_npart; i++){
      vx(i) = vx(i) - sumv(0); //sottraggo alla velocità delle particelle la deriva del sistema -- è la velocità rispetto al riferimento del sistema di particelle
      vy(i) = vy(i) - sumv(1);
      vz(i) = vz(i) - sumv(2);
      sumv2 += vx(i) * vx(i) + vy(i) * vy(i) + vz(i) * vz(i);
    }
    sumv2 /= double(_npart);
    scalef = sqrt(3.0 * _temp / sumv2);   // velocity scale factor
    for (int i=0; i<_npart; i++){
      _particle(i).setvelocity(0, vx(i)*scalef); //alle particelle setto la velocità come quella calcolata rispetto al riferimento delle particelle riscalata per il fattore
      _particle(i).setvelocity(1, vy(i)*scalef);
      _particle(i).setvelocity(2, vz(i)*scalef);
    }
    for (int i=0; i<_npart; i++){
      xold = this->pbc( _particle(i).getposition(0,true) - _particle(i).getvelocity(0)*_delta, 0); //setta le prime posizioni vecchie del sistema eseguendo una prima evoluzione
      yold = this->pbc( _particle(i).getposition(1,true) - _particle(i).getvelocity(1)*_delta, 1);
      zold = this->pbc( _particle(i).getposition(2,true) - _particle(i).getvelocity(2)*_delta, 2);
      _particle(i).setpositold(0, xold);
      _particle(i).setpositold(1, yold);
      _particle(i).setpositold(2, zold);
    }
  }
  return;
}
 
void System :: initialize_properties(){ // Initialize data members used for measurement of properties
 
  string property;
  int index_property = 0; //per ogni proprietà che necessito di misurare aumento di uno l'indice, creando per ogni
  //proprietà un indice che è pari all'indice temporaneo di index_property. L'indice locale serve per spostarsi all'interno
  //del vettore delle misure
  _nprop = 0; //è la lunghezza del vettore delle misure
 
  _measure_penergy  = false; //Defining which properties will be measured
  _measure_kenergy  = false;
  _measure_tenergy  = false;
  _measure_pressure = false;
  _measure_gofr     = false;
  _measure_magnet   = false;
  _measure_cv       = false;
  _measure_chi      = false;
  _measure_pofv     = false;
 
  ifstream input("../INPUT/properties.dat"); //which properties need to be measured?
  if (input.is_open()){
    while ( !input.eof() ){
      input >> property;
      if( property == "POTENTIAL_ENERGY" ){
        ofstream coutp("../OUTPUT/potential_energy.dat");
        coutp << "#     BLOCK:  ACTUAL_PE:     PE_AVE:      ERROR:" << endl; //creo file delle misure
        coutp.close();
        _nprop++; //incremento indice del numero delle proprietà da misurare
        _index_penergy = index_property; //inizializzo posizione della misura dell'energia potenziale nel vettore delle misure
        _measure_penergy = true; //aggiorno flag della misura
        index_property++; //aggiorno indice globale delle proprietà
        _vtail = 0.0; // TO BE FIXED IN EXERCISE 7
      } else if( property == "KINETIC_ENERGY" ){
        ofstream coutk("../OUTPUT/kinetic_energy.dat");
        coutk << "#     BLOCK:   ACTUAL_KE:    KE_AVE:      ERROR:" << endl; //creo file delle misure
        coutk.close();
        _nprop++; //incremento indice del numero dell proprietà da misurare
        _measure_kenergy = true; //aggiorno flag della misura
        _index_kenergy = index_property; //fisso posizione della misura dell'energia cinetica nel vettore delle misure
        index_property++; //incremento indice globale delle proprietà
      } else if( property == "TOTAL_ENERGY" ){
        ofstream coutt("../OUTPUT/total_energy.dat");
        coutt << "#     BLOCK:   ACTUAL_TE:    TE_AVE:      ERROR:" << endl; //creo file misure dell'energia totale
        coutt.close();
        _nprop++;//incremento indice del numero dell proprietà da misurare
        _measure_tenergy = true;//aggiorno flag della misura
        _index_tenergy = index_property;//fisso posizione della misura dell'energia totale nel vettore delle misure
        index_property++; //incremento indice globale delle proprietà
      } else if( property == "TEMPERATURE" ){
        ofstream coutte("../OUTPUT/temperature.dat");
        coutte << "#     BLOCK:   ACTUAL_T:     T_AVE:       ERROR:" << endl; //creo file misura delle temperature
        coutte.close();
        _nprop++;//incremento indice del numero dell proprietà da misurare
        _measure_temp = true;//aggiorno flag della misura
        _index_temp = index_property;//fisso posizione della misura della temperatura nel vettore delle misure
        index_property++;//incremento indice globale delle proprietà
      } else if( property == "PRESSURE" ){
        ofstream coutpr("../OUTPUT/pressure.dat");
        coutpr << "#     BLOCK:   ACTUAL_P:     P_AVE:       ERROR:" << endl; //creo file misure della pressione
        coutpr.close();
        _nprop++;//incremento indice del numero dell proprietà da misurare
        _measure_pressure = true;//aggiorno flag della misura
        _index_pressure = index_property;//fisso posizione della misura della pressione nel vettore delle misure
        index_property++;//incremento indice globale delle proprietà
        _ptail = 0.0; // TO BE FIXED IN EXERCISE 7
      } else if( property == "GOFR" ){
        ofstream coutgr("../OUTPUT/gofr.dat");
        coutgr << "# DISTANCE:     AVE_GOFR:        ERROR:" << endl; //creo file per distanza radiale
        coutgr.close();
        input>>_n_bins; //se vedi nel file properties c'è il numero di bin subito dopo la proprietà da misurare
        _nprop+=_n_bins; //aumenta il numero di proprietà in base alla quantità di bins. Ogni bin è considerata come una proprietà in più
        _bin_size = (_halfside.min() )/(double)_n_bins; //calcola la larghezza del bin --> prendo il lato minore del contenitore e lo dividio per il numero di in
        _measure_gofr = true; //aggiorno il flag della misura
        _index_gofr = index_property;//fisso la posizione nel vettore delle misure. Nota che ho come se avessi un numero di posizioni nele vettore pari al numero di bin
        index_property+= _n_bins;//incremento indice globale delle proprietà considerando il numero di bin
      } else if( property == "POFV" ){
        if(_sim_type > 0){
          cerr << "PROBLEM: DOES NOT MAKE SENSE COMPUTING POFV FOR MC" << endl;
          exit(EXIT_FAILURE);
        }
        // Creo il file di output per la distribuzione delle velocità con intestazione
        ofstream coutpv("../OUTPUT/pofv.dat");
        coutpv << "#  VELOCITY:   AVE_POFV:      ERROR:" << endl; //creo file della distribuzione delle velocità
        coutpv.close();
        input>>_n_bins_v;  // Leggo da properties.dat il numero di bin per POFV
        _nprop += _n_bins_v; /// Incremento il numero totale di proprietà da misurare
        // Imposto la dimensione di ciascun bin come (4·T*)/N_bins_v
        _bin_size_v = 4.0*_temp/(double)_n_bins_v;
 
        _measure_pofv = true;//aggiorno flag della misura
        _index_pofv = index_property;  // Salvo l'indice in _measurement da cui iniziano i bins POFV
        index_property += _n_bins_v; / // Sposto l'indice complessivo considerando tutti i bins POFV
      } else if( property == "MAGNETIZATION" ){
        ofstream coutpr("../OUTPUT/magnetization.dat");
        coutpr << "#     BLOCK:   ACTUAL_M:     M_AVE:       ERROR:" << endl;
        coutpr.close();
        _nprop++;
        double Pressure(int i, int dim); //Calculate pressure
        _measure_magnet = true;
        _index_magnet = index_property;
        index_property++;
      } else if( property == "SPECIFIC_HEAT" ){
        ofstream coutpr("../OUTPUT/specific_heat.dat");
        coutpr << "#     BLOCK:   ACTUAL_CV:    CV_AVE:      ERROR:" << endl;
        coutpr.close();
        _nprop++;
        _measure_cv = true;
        _index_cv = index_property;
        index_property++;
      } else if( property == "SUSCEPTIBILITY" ){
        ofstream coutpr("../OUTPUT/susceptibility.dat");
        coutpr << "#     BLOCK:   ACTUAL_X:     X_AVE:       ERROR:" << endl;
        coutpr.close();
        _nprop++;
        _measure_chi = true;
        _index_chi = index_property;
        index_property++;
      } else if( property == "ENDPROPERTIES" ){
        ofstream coutf;
        coutf.open("../OUTPUT/output.dat",ios::app);
        coutf << "Reading properties completed!" << endl;
        coutf.close();
        break;
      } else cerr << "PROBLEM: unknown property" << endl;
    }
    input.close();
  } else cerr << "PROBLEM: Unable to open properties.dat" << endl;
 
  // according to the number of properties, resize the vectors _measurement,_average,_block_av,_global_av,_global_av2
  _measurement.resize(_nprop); //setto la lunghezza delle misure in base al numero di proprietà calcolate
  _average.resize(_nprop); //setto il vettore delle medie delle proprietà. NOTA le posizioni sono come in measurement
  _block_av.resize(_nprop); // setto il vettore delle medie di blocco
  _global_av.resize(_nprop);
  _global_av2.resize(_nprop);
  _average.zeros(); //setto a zero i tre vettori
  _global_av.zeros();
  _global_av2.zeros();
  _nattempts = 0;
  _naccepted = 0;
  return;
}
 
void System :: finalize(){
  this->write_configuration();
  _rnd.SaveSeed();
  ofstream coutf;
  coutf.open("../OUTPUT/output.dat",ios::app);
  coutf << "Simulation completed!" << endl;
  coutf.close();
  return;
}
 
// Write final configurations as .xyz files in the directory ../OUTPUT/CONFIG/
void System :: write_configuration(){ //salva la configurazione del sistema ad un certo istante
  ofstream coutf;
  if(_sim_type < 2){
    coutf.open("../OUTPUT/CONFIG/config.xyz");
    if(coutf.is_open()){
      coutf << _npart << endl;
      coutf << "#Comment!" << endl;
      for(int i=0; i<_npart; i++){
        coutf << "LJ" << "  "
              << setprecision(17) << _particle(i).getposition(0,true)/_side(0) << "   " // x
              << setprecision(17) << _particle(i).getposition(1,true)/_side(1) << "   " // y
              << setprecision(17) << _particle(i).getposition(2,true)/_side(2) << endl; // z
      }
    } else cerr << "PROBLEM: Unable to open config.xyz" << endl;
    coutf.close();
    coutf.open("../OUTPUT/CONFIG/conf-1.xyz");
    if(coutf.is_open()){
      coutf << _npart << endl;
      coutf << "#Comment!" << endl;
      for(int i=0; i<_npart; i++){
        coutf << "LJ" << "  "
              << setprecision(17) << _particle(i).getposition(0,false)/_side(0) << "   " // x
              << setprecision(17) << _particle(i).getposition(1,false)/_side(1) << "   " // y
              << setprecision(17) << _particle(i).getposition(2,false)/_side(2) << endl; // z
      }
    } else cerr << "PROBLEM: Unable to open conf-1.xyz" << endl;
    coutf.close();
  } else {
    coutf.open("../OUTPUT/CONFIG/config.spin");
    for(int i=0; i<_npart; i++) coutf << _particle(i).getspin() << " ";
    coutf.close();
  }
  return;
}
 
// Write configuration nconf as a .xyz file in directory ../OUTPUT/CONFIG/
void System :: write_XYZ(int nconf){
  ofstream coutf;
  coutf.open("../OUTPUT/CONFIG/config_" + to_string(nconf) + ".xyz");
  if(coutf.is_open()){
    coutf << _npart << endl;
    coutf << "#Comment!" << endl;
    for(int i=0; i<_npart; i++){
      coutf << "LJ" << "  "
            << setw(16) << _particle(i).getposition(0,true)          // x
            << setw(16) << _particle(i).getposition(1,true)          // y
            << setw(16) << _particle(i).getposition(2,true) << endl; // z
    }
  } else cerr << "PROBLEM: Unable to open config.xyz" << endl;
  coutf.close();
  return;
}
 
// Read configuration from a .xyz file in directory ../OUTPUT/CONFIG/
void System :: read_configuration(){
  ifstream cinf;
  cinf.open("../INPUT/CONFIG/config.xyz");
  if(cinf.is_open()){
    string comment;
    string particle;
    double x, y, z;
    int ncoord;
    cinf >> ncoord;
    if (ncoord != _npart){
      cerr << "PROBLEM: conflicting number of coordinates in input.dat & config.xyz not match!" << endl;
      exit(EXIT_FAILURE);
    }
    cinf >> comment;
    for(int i=0; i<_npart; i++){
      cinf >> particle >> x >> y >> z; // units of coordinates in conf.xyz is _side
      _particle(i).setposition(0, this->pbc(_side(0)*x, 0));
      _particle(i).setposition(1, this->pbc(_side(1)*y, 1));
      _particle(i).setposition(2, this->pbc(_side(2)*z, 2));
      _particle(i).acceptmove(); // _x_old = _x_new
    }
  } else cerr << "PROBLEM: Unable to open INPUT file config.xyz"<< endl;
  cinf.close();
  if(_restart and _sim_type > 1){
    int spin;
    cinf.open("../INPUT/CONFIG/config.spin");
    for(int i=0; i<_npart; i++){
      cinf >> spin;
      _particle(i).setspin(spin);
    }
    cinf.close();
  }
  return;
}
 
void System :: block_reset(int blk){ // Reset block accumulators to zero
  ofstream coutf;
  if(blk>0){
    coutf.open("../OUTPUT/output.dat",ios::app);
    coutf << "Block completed: " << blk << endl;
    coutf.close();
  }
  _block_av.zeros();
  return;
}
 
void System :: measure(){ // Measure properties
  _measurement.zeros();
  // POTENTIAL ENERGY, VIRIAL, GOFR ///////////////////////////////////////////
  int bin;
  vec distance;
  distance.resize(_ndim);
  double penergy_temp=0.0, dr; // temporary accumulator for potential energy
  double kenergy_temp=0.0; // temporary accumulator for kinetic energy
  double tenergy_temp=0.0;
  double magnetization=0.0;
  double virial=0.0;
  if (_measure_penergy or _measure_pressure or _measure_gofr) {
    for (int i=0; i<_npart-1; i++){
      for (int j=i+1; j<_npart; j++){
        //calcolo distanza tra coppie di particelle nelle tre dimensioni
        distance(0) = this->pbc( _particle(i).getposition(0,true) - _particle(j).getposition(0,true), 0);
        distance(1) = this->pbc( _particle(i).getposition(1,true) - _particle(j).getposition(1,true), 1);
        distance(2) = this->pbc( _particle(i).getposition(2,true) - _particle(j).getposition(2,true), 2);
        dr = sqrt( dot(distance,distance) ); //modulo della distanza tra particelle
        // GOFR ... TO BE FIXED IN EXERCISE 7
        if(dr < _r_cut){ //faccio il conto fino a che la distanza tra particelle è minore del cutoff del potenziale
          if(_measure_penergy)  penergy_temp += 1.0/pow(dr,12) - 1.0/pow(dr,6); // POTENTIAL ENERGY
          if(_measure_pressure) virial       += 1.0/pow(dr,12) - 0.5/pow(dr,6); // PRESSURE
        }
      }
    }
  }
  if(_measure_pofv){
    
    
    particelle=0;// Contatore di particelle considerate (con v < v_max)
    for(int i = 0; i < _npart; i++){ 
       // Prendo il modulo della velocità della particella i
      if(_particle(i).getvelocity_module() < 4.0*_temp){ //controllo che sia minore di vmax
        particelle ++;
         // Calcolo a quale bin appartiene la velocità: floor(v / bin_size_v)
        int index_occurrence_v = floor((_particle(i).getvelocity_module())/_bin_size_v); 
         // Offsetto con _index_pofv per posizionarmi nella sezione corretta di _measurement
        index_occurrence_v = _index_pofv + index_occurrence_v; 
        // Incremento il conteggio "grezzo" del bin corrispondente
        _measurement(index_occurrence_v)++; 
      }
    }
    if(particelle>0){
       // Normalizzo ciascun bin: dividendo per il numero di particelle totali considerate
      for(int i=0;i<_n_bins_v;i++){
        // Prima normalizzo per il numero di particelle
        _measurement(_index_pofv + i) /= (double)particelle;
        // Poi normalizzo per la larghezza del bin, ottenendo p(v) densità di probabilità
        _measurement(_index_pofv + i) /= _bin_size_v;
      }
    }
  }
  // POTENTIAL ENERGY //////////////////////////////////////////////////////////
  if (_measure_penergy){
    penergy_temp = _vtail + 4.0 * penergy_temp / double(_npart);
    _measurement(_index_penergy) = penergy_temp;
  }
  // KINETIC ENERGY ////////////////////////////////////////////////////////////
  if (_measure_kenergy){ //misura dell'energia cinetica media
    for (int i=0; i<_npart; i++) kenergy_temp += 0.5 * dot( _particle(i).getvelocity() , _particle(i).getvelocity() );
    kenergy_temp /= double(_npart);
    _measurement(_index_kenergy) = kenergy_temp; //il vettore measurment viene costantemente aggiornato in quella posizione per ogni misura istantanea
  }
  // TOTAL ENERGY (kinetic+potential) //////////////////////////////////////////
  if (_measure_tenergy){
    if (_sim_type < 2) _measurement(_index_tenergy) = kenergy_temp + penergy_temp;
    else {
      double s_i, s_j;
      for (int i=0; i<_npart; i++){
        s_i = double(_particle(i).getspin());
        s_j = double(_particle(this->pbc(i+1)).getspin());
        tenergy_temp += - _J * s_i * s_j - 0.5 * _H * (s_i + s_j);
      }
      tenergy_temp /= double(_npart);
      _measurement(_index_tenergy) = tenergy_temp;
    }
  }
  // TEMPERATURE ///////////////////////////////////////////////////////////////
  if (_measure_temp and _measure_kenergy) _measurement(_index_temp) = (2.0/3.0) * kenergy_temp;
  // PRESSURE //////////////////////////////////////////////////////////////////
  if (_measure_pressure) _measurement[_index_pressure] = _rho * (2.0/3.0) * kenergy_temp + (_ptail*_npart + 48.0*virial/3.0)/_volume;
  // MAGNETIZATION /////////////////////////////////////////////////////////////
// TO BE FIXED IN EXERCISE 6
  // SPECIFIC HEAT /////////////////////////////////////////////////////////////
// TO BE FIXED IN EXERCISE 6
  // SUSCEPTIBILITY ////////////////////////////////////////////////////////////
// TO BE FIXED IN EXERCISE 6
 
  _block_av += _measurement; //Update block accumulators -- E' un vettore le cui posizioni sono uguali a quelle di meausremente e viene ad ogni misura sommato
 
  return;
}
 
void System :: averages(int blk){ //come input metto il numero del blk-esimo blocco
 
  ofstream coutf;
  double average, sum_average, sum_ave2;
  vec bins_average, sum_bins_average, sum2_bins_average;
 
  _average     = _block_av / double(_nsteps); //E' un vettore di medie di misure del blocco corrente
  _global_av  += _average; //viene sommato il vettore delle medie nel blocco corrente
  _global_av2 += _average % _average; // % -> element-wise multiplication -- vettore per calcolare le incertezze
 
  // POTENTIAL ENERGY //////////////////////////////////////////////////////////
  if (_measure_penergy){
    coutf.open("../OUTPUT/potential_energy.dat",ios::app);
    average  = _average(_index_penergy);
    sum_average = _global_av(_index_penergy);
    sum_ave2 = _global_av2(_index_penergy);
    coutf << setw(12) << blk
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  // KINETIC ENERGY ////////////////////////////////////////////////////////////
  if (_measure_kenergy){
    coutf.open("../OUTPUT/kinetic_energy.dat",ios::app); //questo file è stato creato quando si è inizialiaato il sistema
    average  = _average(_index_kenergy); //selezione la media all'indice dell'energia cinetica
    sum_average = _global_av(_index_kenergy); //è la somma delle medie -- seleziono quella alla posizione dell'energia cinetica
    sum_ave2 = _global_av2(_index_kenergy);//è la somma dei quadrati delle medie -- seleziono quella alla posizione dell'energia cinetica
    coutf << setw(12) << blk //ricorda setw() setta la lunghezza minima dell'output --> verrano forzati 12 caratteri in scrittura
          << setw(12) << average //stampa la media del blocco numero blk
          << setw(12) << sum_average/double(blk) //stampa la media delle medie fino al blocco blk
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl; //calcola l'errore sul blocco blk
    coutf.close();
  }
  // TOTAL ENERGY //////////////////////////////////////////////////////////////
  if (_measure_tenergy){
    coutf.open("../OUTPUT/total_energy.dat",ios::app);
    average  = _average(_index_tenergy);
    sum_average = _global_av(_index_tenergy);
    sum_ave2 = _global_av2(_index_tenergy);
    coutf << setw(12) << blk
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  // TEMPERATURE ///////////////////////////////////////////////////////////////
  if (_measure_temp){
    coutf.open("../OUTPUT/temperature.dat",ios::app);
    average  = _average(_index_temp);
    sum_average = _global_av(_index_temp);
    sum_ave2 = _global_av2(_index_temp);
    coutf << setw(12) << blk
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  // PRESSURE //////////////////////////////////////////////////////////////////
  if (_measure_pressure){
    coutf.open("../OUTPUT/pressure.dat",ios::app);
    average  = _average(_index_pressure); //media di blocco
    sum_average = _global_av(_index_pressure); //media delle medie al blocco attuale
    sum_ave2 = _global_av2(_index_pressure);//calcolo errore -- metto somma delle medie fino al blocco corrente,
    // somma delle medie quadre fino al blocco corrente e il blocco attuale
    coutf << setw(12) << blk
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  // GOFR //////////////////////////////////////////////////////////////////////
  // TO BE FIXED IN EXERCISE 7
  // POFV //////////////////////////////////////////////////////////////////////
  if(_measure_pofv && blk == _nblocks){ // stampa solo all'ultimo blocco
     // Stampo i risultati di POFV solo all’ultimo blocco
    ofstream coutf("../OUTPUT/pofv.dat", ios::app);
    for(int i = 0; i < _n_bins_v; i++){
       // Centro del bin: (i + 0.5) * bin_size_v
      double bin_center = (i + 0.5) * _bin_size_v;
       // Calcolo la media progressiva: media delle medie blocco su tutti i blocchi
      double avg = _global_av(_index_pofv + i) / double(blk);
      // Calcolo l’errore statistico dalla somma delle medie al quadrato
      double err = this->error(_global_av(_index_pofv + i), _global_av2(_index_pofv + i), blk);
      
      coutf << setw(12) << bin_center
            << setw(12) << avg
            << setw(12) << err << endl;
    }
    coutf.close();
  }
  
  // MAGNETIZATION /////////////////////////////////////////////////////////////
  // TO BE FIXED IN EXERCISE 6
  // SPECIFIC HEAT /////////////////////////////////////////////////////////////
  // TO BE FIXED IN EXERCISE 6
  // SUSCEPTIBILITY ////////////////////////////////////////////////////////////
  // TO BE FIXED IN EXERCISE 6
  // ACCEPTANCE ////////////////////////////////////////////////////////////////
  double fraction;
  coutf.open("../OUTPUT/acceptance.dat",ios::app);
  if(_nattempts > 0) fraction = double(_naccepted)/double(_nattempts);
  else fraction = 0.0;
  coutf << setw(12) << blk << setw(12) << fraction << endl;
  coutf.close();
  
  return;
}
 
double System :: error(double acc, double acc2, int blk){ //acc è la somma delle medie fino al blocco blk, acc2 è la sommma delle medie dei blocchi quadri fino al blocco blk
  if(blk <= 1) return 0.0; //se il blocco è il primo l'incertezza la metto a zero
  else return sqrt( fabs(acc2/double(blk) - pow( acc/double(blk) ,2) )/double(blk) ); //altrimenti calcolo la deviazione standard sulla media dei blocchi
}
 
int System :: get_nbl(){ //ritorna numero dei blocchi
  return _nblocks;
}
 
int System :: get_nsteps(){//ritorna numero di step per blocco
  return _nsteps;
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
 
 