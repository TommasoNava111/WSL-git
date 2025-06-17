#include <cmath>
#include <cstdlib>
#include <string>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "system.h"

using namespace std;
using namespace arma;

void System::step(){ // Perform a simulation step
  if(_sim_type == 0) 
    this->Verlet();  // Perform a MD step
  else 
    for(int i=0; i<_npart; i++) 
      this->move(int(_rnd.Rannyu()*_npart)); // Perform a MC step on a randomly chosen particle
  _nattempts += _npart; // update number of attempts performed on the system
  return;
}

void System::Verlet(){
  double xnew, ynew, znew;
  // Force acting on each particle
  for(int i=0; i<_npart; i++){
    _fx(i) = this->Force(i,0);
    _fy(i) = this->Force(i,1);
    _fz(i) = this->Force(i,2);
  }
  // Verlet integration scheme for each particle
  for(int i=0; i<_npart; i++){
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

double System::Force(int i, int dim){
  double f=0.0, dr;
  vec distance;
  distance.resize(_ndim);
  for (int j=0; j<_npart; j++){
    if(i != j){
      distance(0) = this->pbc( _particle(i).getposition(0,true) - _particle(j).getposition(0,true), 0);
      distance(1) = this->pbc( _particle(i).getposition(1,true) - _particle(j).getposition(1,true), 1);
      distance(2) = this->pbc( _particle(i).getposition(2,true) - _particle(j).getposition(2,true), 2);
      dr = sqrt(dot(distance,distance));
      if(dr < _r_cut){
        f += distance(dim) * (48.0/pow(dr,14) - 24.0/pow(dr,8));
      }
    }
  }
  return f;
}

void System::move(int i){ // Provo a "muovere" (cioè a cambiare lo spin) del sito i
  if(_sim_type == 3){ // Gibbs sampler per Ising (sim_type==3)
    // Calcolo lo spin dei vicini (i-1 e i+1) con condizioni di periodicità:
    int s_left  = _particle(this->pbc(i-1)).getspin();
    int s_right = _particle(this->pbc(i+1)).getspin();
    // Probabilità condizionata di avere spin +1:
    double p_up = 1.0/(1.0 + exp(-2.0*_beta*(_J*(s_left+s_right) + _H)));
    // Genero un numero casuale in [0,1). Se < p_up, imposto spin = +1, altrimenti -1.
    if(_rnd.Rannyu() < p_up) 
      _particle(i).setspin(1);
    else 
      _particle(i).setspin(-1);
    _naccepted++; // La mossa è sempre "accettata" nel Gibbs sampler
  } else { // Metropolis per Ising (sim_type==2) o altri sistemi (sim_type==1)
    if(_sim_type == 1){       // LJ system
      vec shift(_ndim);       // Memorizza la proposta di traslazione
      for(int j=0; j<_ndim; j++){
        shift(j) = _rnd.Rannyu(-1.0,1.0) * _delta; // distribuzione uniforme in [-_delta, _delta)
      }
      _particle(i).translate(shift, _side);  // Chiamata a Particle::translate
      if(this->metro(i)){ // Valutazione del criterio di Metropolis
        _particle(i).acceptmove();
        _naccepted++;
      } else 
        _particle(i).moveback(); // Se rifiutato, ripristina la configurazione precedente
    } else {                  // Ising 1D con Metropolis (sim_type==2)
      if(this->metro(i)){     // Valutazione del criterio di Metropolis per un flip di spin
        _particle(i).flip();  // Se accettato, viene effettuato il flip
        _naccepted++;
      } // Se metro(i) restituisce false, non faccio nulla: lo spin resta com’era
    }
  }
  return;
}

bool System::metro(int i){ // Algoritmo di Metropolis
  bool decision = false;
  double delta_E, acceptance;
  if(_sim_type == 1) 
    delta_E = this->Boltzmann(i,true) - this->Boltzmann(i,false);
  else 
   // Calcolo della variazione di energia ΔE se flippassi lo spin i:
    delta_E = 2.0 * _particle(i).getspin() * 
              ( _J * (_particle(this->pbc(i-1)).getspin() + _particle(this->pbc(i+1)).getspin() ) + _H );
  // Accettazione di Metropolis: exp(-β ΔE)
  acceptance = exp(-_beta*delta_E);
  // Confronto con un numero casuale [0,1):
  if(_rnd.Rannyu() < acceptance ) 
    decision = true; // Step di accettazione di Metropolis
  return decision;
}

double System::Boltzmann(int i, bool xnew){
  double energy_i=0.0;
  double dx, dy, dz, dr;
  for (int j=0; j<_npart; j++){
    if(j != i){
      dx = this->pbc(_particle(i).getposition(0,xnew) - _particle(j).getposition(0,1), 0);
      dy = this->pbc(_particle(i).getposition(1,xnew) - _particle(j).getposition(1,1), 1);
      dz = this->pbc(_particle(i).getposition(2,xnew) - _particle(j).getposition(2,1), 2);
      dr = sqrt(dx*dx + dy*dy + dz*dz);
      if(dr < _r_cut){
        energy_i += 1.0/pow(dr,12) - 1.0/pow(dr,6);
      }
    }
  }
  return 4.0 * energy_i;
}

double System::pbc(double position, int i){ // Applica le condizioni periodiche
  return position - _side(i) * rint(position / _side(i));
}

int System::pbc(int i){ // Versione per gli spin (in 1D)
  if(i >= _npart) i = i - _npart;
  else if(i < 0)  i = i + _npart;
  return i;
} 

void System::initialize(){ // Inizializza l'oggetto System leggendo i file in ../INPUT/
  int p1, p2; // Legge dal file ../INPUT/Primes
  ifstream Primes("../INPUT/Primes");
  Primes >> p1 >> p2;
  Primes.close();
  int seed[4]; // Legge il seed per il RNG
  ifstream Seed("../INPUT/seed.in");
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  _rnd.SetRandom(seed, p1, p2);

  ofstream couta("../OUTPUT/acceptance.dat"); // Inizializza l'header per acceptance.dat
  couta << "#   N_BLOCK:  ACCEPTANCE:" << endl;
  couta.close();

  ifstream input("../INPUT/input.dat"); // Inizia la lettura di input.dat
  ofstream coutf;
  coutf.open("../OUTPUT/output.dat");
  string property;
  double delta;
  while ( !input.eof() ){
    input >> property;
    if( property == "SIMULATION_TYPE" ){
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
      else if(_sim_type == 3) coutf << "ISING 1D MONTE CARLO (GIBBS) SIMULATION"   << endl;
    } else if( property == "RESTART" ){
      input >> _restart;
    } else if( property == "TEMP" ){
      input >> _temp;
      _beta = 1.0/_temp;
      coutf << "TEMPERATURE= " << _temp << endl;
    } else if( property == "NPART" ){
      input >> _npart;
      _fx.resize(_npart);
      _fy.resize(_npart);
      _fz.resize(_npart);
      _particle.set_size(_npart);
      for(int i=0; i<_npart; i++){ 
        _particle(i).initialize();
        if(_rnd.Rannyu() > 0.5) 
          _particle(i).flip(); // Randomizza la configurazione degli spin
      }
      coutf << "NPART= " << _npart << endl;
    } else if( property == "RHO" ){
      input >> _rho;
      _volume = _npart/_rho;
      _side.resize(_ndim);
      _halfside.resize(_ndim);
      double side = pow(_volume, 1.0/3.0);
      for(int i=0; i<_ndim; i++) _side(i) = side;
      _halfside = 0.5 * _side;
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
    } else 
      cerr << "PROBLEM: unknown input" << endl;
  }
  input.close();
  this->read_configuration();
  if(_sim_type==0) 
    this->initialize_velocities();
  coutf << "System initialized!" << endl;
  coutf.close();
  return;
}

void System::initialize_velocities(){
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
        cinf >> particle >> xold >> yold >> zold; // coordinate normalizzate in unità di _side
        _particle(i).setpositold(0, this->pbc(_side(0)*xold, 0));
        _particle(i).setpositold(1, this->pbc(_side(1)*yold, 1));
        _particle(i).setpositold(2, this->pbc(_side(2)*zold, 2));
      }
    } else 
      cerr << "PROBLEM: Unable to open INPUT file conf-1.xyz"<< endl;
    cinf.close();
  } else {
    vec vx(_npart), vy(_npart), vz(_npart);
    vec sumv(_ndim);
    sumv.zeros();
    for (int i=0; i<_npart; i++){
      vx(i) = _rnd.Gauss(0.,sqrt(_temp));
      vy(i) = _rnd.Gauss(0.,sqrt(_temp));
      vz(i) = _rnd.Gauss(0.,sqrt(_temp));
      sumv(0) += vx(i);
      sumv(1) += vy(i);
      sumv(2) += vz(i);
    }
    for (int idim=0; idim<_ndim; idim++) 
      sumv(idim) = sumv(idim)/double(_npart);
    double sumv2 = 0.0, scalef;
    for (int i=0; i<_npart; i++){
      vx(i) = vx(i) - sumv(0);
      vy(i) = vy(i) - sumv(1);
      vz(i) = vz(i) - sumv(2);
      sumv2 += vx(i) * vx(i) + vy(i) * vy(i) + vz(i) * vz(i);
    }
    sumv2 /= double(_npart);
    scalef = sqrt(3.0 * _temp / sumv2);   // fattore di scala per le velocità 
    for (int i=0; i<_npart; i++){
      _particle(i).setvelocity(0, vx(i)*scalef);
      _particle(i).setvelocity(1, vy(i)*scalef);
      _particle(i).setvelocity(2, vz(i)*scalef);
    }
    for (int i=0; i<_npart; i++){
      xold = this->pbc( _particle(i).getposition(0,true) - _particle(i).getvelocity(0)*_delta, 0);
      yold = this->pbc( _particle(i).getposition(1,true) - _particle(i).getvelocity(1)*_delta, 1);
      zold = this->pbc( _particle(i).getposition(2,true) - _particle(i).getvelocity(2)*_delta, 2);
      _particle(i).setpositold(0, xold);
      _particle(i).setpositold(1, yold);
      _particle(i).setpositold(2, zold);
    }
  }
  return;
}

void System::initialize_properties(){ // Inizializza le grandezze per le misure
  string property;
  int index_property = 0;
  _nprop = 0;

  _measure_penergy  = false; // Definisce quali proprietà verranno misurate
  _measure_kenergy  = false;
  _measure_tenergy  = false;
  _measure_pressure = false;
  _measure_gofr     = false;
  _measure_magnet   = false;
  _measure_cv       = false;
  _measure_chi      = false;

  // Se la simulazione è Ising (_sim_type>=2) definiamo il suffisso per i file
  string suffix = "";
  if(_sim_type >= 2){
    suffix = (_sim_type == 2) ? "_metropolis" : "_gibbs";
  }

  ifstream input("../INPUT/properties.dat");
  if (input.is_open()){
    while ( !input.eof() ){
      input >> property;
      if( property == "POTENTIAL_ENERGY" ){
        ofstream coutp("../OUTPUT/potential_energy.dat");
        coutp << "#     BLOCK:  ACTUAL_PE:     PE_AVE:      ERROR:" << endl;
        coutp.close();
        _nprop++;
        _index_penergy = index_property;
        _measure_penergy = true;
        index_property++;
        _vtail = 0.0; // TO BE FIXED IN EXERCISE 7
      } else if( property == "KINETIC_ENERGY" ){
        ofstream coutk("../OUTPUT/kinetic_energy.dat");
        coutk << "#     BLOCK:   ACTUAL_KE:    KE_AVE:      ERROR:" << endl;
        coutk.close();
        _nprop++;
        _measure_kenergy = true;
        _index_kenergy = index_property;
        index_property++;
      } else if( property == "TOTAL_ENERGY" ){
        
        _nprop++;
        _measure_tenergy = true;
        _index_tenergy = index_property;
        index_property++;
      } else if( property == "TEMPERATURE" ){
        ofstream coutte("../OUTPUT/temperature.dat");
        coutte << "#     BLOCK:   ACTUAL_T:     T_AVE:       ERROR:" << endl;
        coutte.close();
        _nprop++;
        _measure_temp = true;
        _index_temp = index_property;
        index_property++;
      } else if( property == "PRESSURE" ){
        ofstream coutpr("../OUTPUT/pressure.dat");
        coutpr << "#     BLOCK:   ACTUAL_P:     P_AVE:       ERROR:" << endl;
        coutpr.close();
        _nprop++;
        _measure_pressure = true;
        _index_pressure = index_property;
        index_property++;
        _ptail = 0.0; // TO BE FIXED IN EXERCISE 7
      } else if( property == "GOFR" ){
        ofstream coutgr("../OUTPUT/gofr.dat");
        coutgr << "# DISTANCE:     AVE_GOFR:        ERROR:" << endl;
        coutgr.close();
        input >> _n_bins;
        _nprop += _n_bins;
        _bin_size = (_halfside.min())/(double)_n_bins;
        _measure_gofr = true;
        _index_gofr = index_property;
        index_property += _n_bins;
      } else if( property == "POFV" ){
        if(_sim_type > 0){
          cerr << "PROBLEM: DOES NOT MAKE SENSE COMPUTING POFV FOR MC" << endl;
          exit(EXIT_FAILURE);
        }
        ofstream coutpv("../OUTPUT/pofv.dat");
        coutpv << "# VELOCITY:     AVE_POFV:        ERROR:" << endl;
        coutpv.close();
        input >> _n_bins_v;
        _nprop += _n_bins_v;
        _bin_size_v = 4.0*_temp/(double)_n_bins_v; // TO BE FIXED IN EXERCISE 4
        _measure_pofv = true;
        _index_pofv = index_property;
        index_property += _n_bins_v;
      } else if( property == "MAGNETIZATION" ){
        
        _nprop++;
        _measure_magnet = true;
        _index_magnet = index_property;
        index_property++;
      } else if( property == "SPECIFIC_HEAT" ){
        
        _nprop++;
        _measure_cv = true;
        _index_cv = index_property;
        index_property++;
      } else if( property == "SUSCEPTIBILITY" ){
        
         _nprop++;
        _measure_chi = true;
        _index_chi = index_property;
        index_property++;
      } else if( property == "ENDPROPERTIES" ){
        ofstream coutf;
        coutf.open("../OUTPUT/output.dat", ios::app);
        coutf << "Reading properties completed!" << endl;
        coutf.close();
        break;
      } else 
        cerr << "PROBLEM: unknown property" << endl;
    }
    input.close();
  } else 
    cerr << "PROBLEM: Unable to open properties.dat" << endl;

  // Ridimensiona i vettori per le misurazioni in base al numero di proprietà
  _measurement.resize(_nprop);
  _average.resize(_nprop);
  _block_av.resize(_nprop);
  _global_av.resize(_nprop);
  _global_av2.resize(_nprop);
  _average.zeros();
  _global_av.zeros();
  _global_av2.zeros();
  _nattempts = 0;
  _naccepted = 0;
  return;
}

void System::finalize(){
  this->write_configuration();
  _rnd.SaveSeed();
  ofstream coutf;
  coutf.open("../OUTPUT/output.dat", ios::app);
  coutf << "Simulation completed!" << endl;
  coutf.close();
  return;
}

// Scrive la configurazione finale in file .xyz o .spin in ../OUTPUT/CONFIG/
void System::write_configuration(){
  ofstream coutf;
  if(_sim_type < 2){
    coutf.open("../OUTPUT/CONFIG/config.xyz");
    if(coutf.is_open()){
      coutf << _npart << endl;
      coutf << "#Comment!" << endl;
      for(int i=0; i<_npart; i++){
        coutf << "LJ" << "  " 
              << setprecision(17) << _particle(i).getposition(0,true)/_side(0) << "   " 
              << setprecision(17) << _particle(i).getposition(1,true)/_side(1) << "   " 
              << setprecision(17) << _particle(i).getposition(2,true)/_side(2) << endl;
      }
    } else 
      cerr << "PROBLEM: Unable to open config.xyz" << endl;
    coutf.close();
    coutf.open("../OUTPUT/CONFIG/conf-1.xyz");
    if(coutf.is_open()){
      coutf << _npart << endl;
      coutf << "#Comment!" << endl;
      for(int i=0; i<_npart; i++){
        coutf << "LJ" << "  "
              << setprecision(17) << _particle(i).getposition(0,false)/_side(0) << "   " 
              << setprecision(17) << _particle(i).getposition(1,false)/_side(1) << "   " 
              << setprecision(17) << _particle(i).getposition(2,false)/_side(2) << endl;
      }
    } else 
      cerr << "PROBLEM: Unable to open conf-1.xyz" << endl;
    coutf.close();
  } else {
    coutf.open("../OUTPUT/CONFIG/config.spin");
    for(int i=0; i<_npart; i++) 
      coutf << _particle(i).getspin() << " ";
    coutf.close();
  }
  return;
}

// Scrive la configurazione nconf in un file .xyz in ../OUTPUT/CONFIG/
void System::write_XYZ(int nconf){
  ofstream coutf;
  coutf.open("../OUTPUT/CONFIG/config_" + to_string(nconf) + ".xyz");
  if(coutf.is_open()){
    coutf << _npart << endl;
    coutf << "#Comment!" << endl;
    for(int i=0; i<_npart; i++){
      coutf << "LJ" << "  " 
            << setw(16) << _particle(i).getposition(0,true)
            << setw(16) << _particle(i).getposition(1,true)
            << setw(16) << _particle(i).getposition(2,true) << endl;
    }
  } else 
    cerr << "PROBLEM: Unable to open config.xyz" << endl;
  coutf.close();
  return;
}

// Legge la configurazione da un file .xyz in ../INPUT/CONFIG/
void System::read_configuration(){
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
      cinf >> particle >> x >> y >> z; // coordinate normalizzate in unità di _side
      _particle(i).setposition(0, this->pbc(_side(0)*x, 0));
      _particle(i).setposition(1, this->pbc(_side(1)*y, 1));
      _particle(i).setposition(2, this->pbc(_side(2)*z, 2));
      _particle(i).acceptmove(); // _x_old = _x_new
    }
  } else 
    cerr << "PROBLEM: Unable to open INPUT file config.xyz"<< endl;
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

void System::block_reset(int blk){ // Reset dei contatori di blocco
  ofstream coutf;
  if(blk>0){
    coutf.open("../OUTPUT/output.dat", ios::app);
    coutf << "Block completed: " << blk << endl;
    coutf.close();
  }
  _block_av.zeros();
  return;
}

void System::measure(){ // Misura le proprietà
  _measurement.zeros(); // azzero il vettore di questa singola misura
  // Questa funzione viene chiamata dopo ogni step(), e accumula in _measurement[]:
  //  - energia totale per spin (index_tenergy)
  //  - E_total^2, per la capacità termica (index_cv)
  //  - β (Σ s_i)^2 / N, per la suscettività (index_chi)
  //  - (Σ s_i)/N, per la magnetizzazione (index_magnet)
  // Poi aggiunge _measurement a _block_av.
  // Per il sistema LJ: POTENTIAL ENERGY, VIRIAL e GOFR
  int bin;
  vec distance;
  distance.resize(_ndim);
  double penergy_temp=0.0, dr;
  double kenergy_temp=0.0;
  double tenergy_temp=0.0; // accumulatore energia totale (non diviso per spin)
  double magnetization=0.0; // accumulatore somma degli spin
  double virial=0.0;
  
  if (_measure_penergy or _measure_pressure or _measure_gofr) {
    for (int i=0; i<_npart-1; i++){
      for (int j=i+1; j<_npart; j++){
        distance(0) = this->pbc( _particle(i).getposition(0,true) - _particle(j).getposition(0,true), 0);
        distance(1) = this->pbc( _particle(i).getposition(1,true) - _particle(j).getposition(1,true), 1);
        distance(2) = this->pbc( _particle(i).getposition(2,true) - _particle(j).getposition(2,true), 2);
        dr = sqrt(dot(distance,distance));
        if(dr < _r_cut){
          if(_measure_penergy)  
            penergy_temp += 1.0/pow(dr,12) - 1.0/pow(dr,6);
          if(_measure_pressure) 
            virial += 1.0/pow(dr,12) - 0.5/pow(dr,6);
        }
      }
    }
  }
  
  // POTENTIAL ENERGY per LJ
  if (_measure_penergy){
    penergy_temp = _vtail + 4.0 * penergy_temp / double(_npart);
    _measurement(_index_penergy) = penergy_temp;
  }
  // KINETIC ENERGY per LJ
  if (_measure_kenergy){
    for (int i=0; i<_npart; i++) 
      kenergy_temp += 0.5 * dot( _particle(i).getvelocity() , _particle(i).getvelocity() ); 
    kenergy_temp /= double(_npart);
    _measurement(_index_kenergy) = kenergy_temp;
  }
  // Se sto misurando energia totale (PING: _measure_tenergy==true se sim_type>=2)
  if (_measure_tenergy){
    if (_sim_type < 2) 
      _measurement(_index_tenergy) = kenergy_temp + penergy_temp;
    else { 
      double s_i, s_j;
      // Scorro tutti gli spin per calcolare energia di accoppiamento tra i ed i+1:
      for (int i=0; i<_npart; i++){
        s_i = double(_particle(i).getspin());
        s_j = double(_particle(this->pbc(i+1)).getspin());
        // Hamiltoniana 1D per spin:   -J s_i s_{i+1}  - (H/2) (s_i + s_{i+1})
      // Se sommassi su i da 1 a N, i termini -(H/2) (s_i + s_{i+1}) diventano -H s_i (ogni spin appare due volte)
        tenergy_temp += - _J * s_i * s_j - 0.5 * _H * (s_i + s_j);
      }
      // Ora tenergy_temp è la somma totale su tutti gli spin (notare: non è ancora diviso per N)
      tenergy_temp /= double(_npart);
       // Salvo l’energia "per spin" in index_tenergy:
      _measurement(_index_tenergy) = tenergy_temp;
    }
  }
  // TEMPERATURE
  if (_measure_temp and _measure_kenergy) 
    _measurement(_index_temp) = (2.0/3.0) * kenergy_temp;
  // PRESSURE
  if (_measure_pressure) 
    _measurement(_index_pressure) = _rho * (2.0/3.0) * kenergy_temp + (_ptail*_npart + 48.0*virial/3.0)/_volume;

  // --- Misurazioni per il modello Ising (sim_type>=2) ---
  if (_sim_type >= 2) {
    // SPECIFICHE: l'energia interna per spin (già in _index_tenergy)
    // Aggiungo l'energia al quadrato, per la capacità termica:
    if(_measure_cv) {
      double E_total = tenergy_temp * _npart; // Energia totale, non per spin
      _measurement(_index_cv) = E_total * E_total;
    }
    // SUSCEPTIBILITY: beta * (<(somma degli spin)^2>) (normalizzata per spin)
    if(_measure_chi) {
      double sum_s = 0.0; // somma degli spin, usata per suscettività
      for (int i=0; i<_npart; i++){
        sum_s += _particle(i).getspin();
      }
      // Suscettività per spin = β ( Σ s_i )^2 / N
      _measurement(_index_chi) = _beta * (sum_s * sum_s) / _npart;
    }
    // MAGNETIZATION per spin:
    if(_measure_magnet) {
      double magnet = 0.0;
      for (int i=0; i<_npart; i++){
        magnet += _particle(i).getspin();
      }
       // Magnetizzazione per spin = (Σ s_i) / N
      _measurement(_index_magnet) = magnet / double(_npart);
    }
  }
  
  _block_av += _measurement; // Aggiorna accumulatore del blocco
  return;
}

void System::averages(int blk){
  ofstream coutf;
  double average, sum_average, sum_ave2;

  // Definiamo il suffisso per le proprietà Ising se applicabile
  string suffix = "";
  if(_sim_type >= 2){
    suffix = (_sim_type == 2) ? "_metropolis" : "_gibbs";
  }

  _average     = _block_av / double(_nsteps);
  _global_av  += _average;
  _global_av2 += _average % _average; // Element-wise multiplication

  // POTENTIAL ENERGY
  if (_measure_penergy){
    coutf.open("../OUTPUT/potential_energy.dat", ios::app);
    average  = _average(_index_penergy);
    sum_average = _global_av(_index_penergy);
    sum_ave2 = _global_av2(_index_penergy);
    coutf << setw(12) << blk 
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  // KINETIC ENERGY
  if (_measure_kenergy){
    coutf.open("../OUTPUT/kinetic_energy.dat", ios::app);
    average  = _average(_index_kenergy);
    sum_average = _global_av(_index_kenergy);
    sum_ave2 = _global_av2(_index_kenergy);
    coutf << setw(12) << blk
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  // TOTAL ENERGY
  if (_measure_tenergy){
    
    average  = _average(_index_tenergy);
    sum_average = _global_av(_index_tenergy);
    sum_ave2 = _global_av2(_index_tenergy);
   
  }
  // TEMPERATURE
  if (_measure_temp){
    coutf.open("../OUTPUT/temperature.dat", ios::app);
    average  = _average(_index_temp);
    sum_average = _global_av(_index_temp);
    sum_ave2 = _global_av2(_index_temp);
    coutf << setw(12) << blk
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  // PRESSURE
  if (_measure_pressure){
    coutf.open("../OUTPUT/pressure.dat", ios::app);
    average  = _average(_index_pressure);
    sum_average = _global_av(_index_pressure);
    sum_ave2 = _global_av2(_index_pressure);
    coutf << setw(12) << blk
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  //POFV
  if (_measure_pofv && blk == _nblocks) {  // Scrivi solo all'ultimo blocco
    coutf.open("../OUTPUT/pofv" + suffix + ".dat");
    double total_samples = _npart * _nblocks * _nsteps;
    double norm_factor = 1.0 / (total_samples * _bin_size_v );
    for (int i = 0; i < _n_bins_v; ++i) {
        double v_center = (i + 0.5) * _bin_size_v;
        double sum_av = _global_av(_index_pofv + i) * _nsteps;
        double sum_av2 = _global_av2(_index_pofv + i) * _nsteps * _nsteps;
        double avg = sum_av * norm_factor;
        double variance = (sum_av2 / _nblocks - pow(sum_av / _nblocks, 2)) / _nblocks;
        double err = sqrt(variance) * norm_factor;
        coutf << setw(12) << v_center
            << setw(12) << avg
            << setw(12) << err << endl;
    }
    coutf.close();
  }
  // --- Proprietà per il modello Ising ---
  // MAGNETIZATION
  if (_measure_magnet){
    
    average  = _average(_index_magnet);
    sum_average = _global_av(_index_magnet);
    sum_ave2 = _global_av2(_index_magnet);
    
  }
  // SPECIFIC HEAT
  if (_measure_cv) {
    
    // Calcola la varianza dell'energia totale
    double avg_Etot = _global_av(_index_tenergy) / double(blk) * _npart; 
    double avg_Etot2 = _global_av(_index_cv) / double(blk);
    double Cv = (_beta * _beta) * (avg_Etot2 - avg_Etot * avg_Etot) / _npart;
    // Calcola l'errore
    double sum_Etot2 = _global_av(_index_cv);
    double sum_Etot2_sq = _global_av2(_index_cv);
    double err = this->error(sum_Etot2, sum_Etot2_sq, blk) * (_beta * _beta) / _npart;
   
  }
  // SUSCEPTIBILITY
  if (_measure_chi){
   
    average  = _average(_index_chi);
    sum_average = _global_av(_index_chi);
    sum_ave2 = _global_av2(_index_chi);
   
  }
  
  // ACCEPTANCE (per MC, solo per i sistemi non dinamici)
  double fraction;
  coutf.open("../OUTPUT/acceptance.dat", ios::app);
  if(_nattempts > 0) 
    fraction = double(_naccepted)/double(_nattempts);
  else 
    fraction = 0.0; 
  coutf << setw(12) << blk << setw(12) << fraction << endl;
  coutf.close();
  
  return;
}

double System::error(double acc, double acc2, int blk){
  if(blk <= 1) 
    return 0.0;
  else 
    return sqrt(fabs(acc2/double(blk) - pow(acc/double(blk),2)) / double(blk));
}

int System::get_nbl(){
  return _nblocks;
}

int System::get_nsteps(){
  return _nsteps;
}
void System::set_T(double T) {
  _temp = T;
  _beta = 1.0 / T;
}

void System::set_H(double H) {
  _H = H;
}

void System::reset_measurements() {
  _global_av.zeros();
  _global_av2.zeros();
  _block_av.zeros();
  _naccepted = 0;
  _nattempts = 0;
}
