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
 * particle.cpp: Implementazione della classe Particle          *
 *                                                              *
 * Questo file definisce i metodi associati agli oggetti        *
 * Particle, che rappresentano le singole particelle in         *
 * simulazioni di dinamica molecolare (Lennard-Jones) o nel      *
 * modello di Ising 1D.                                          *
 *                                                              *
 * Struttura del file:                                           *
 *                                                              *
 * 1) Inclusioni                                                *
 *    - <iostream>: per eventuali operazioni di I/O di debug.   *
 *    - <math.h>: per funzioni matematiche come sqrt, pow.      *
 *    - "particle.h": dichiarazione della classe Particle.      *
 *                                                              *
 * 2) void Particle::initialize()                               *
 *    - Inizializza lo spin a +1 (utile per Ising).             *
 *    - Alloca e ridimensiona i vettori privati:                 *
 *        • _x, _xold: dimensione _ndim (coordinate corrente e    *
 *          “vecchie”).                                          *
 *        • _v: vettore di velocità (dimensione _ndim).           *
 *    - Non ritorna nulla (void).                                *
 *                                                              *
 * 3) void Particle::translate(vec delta, vec side)             *
 *    - Applica uno spostamento “delta” alle coordinate _x.      *
 *      Per ogni dimensione i (0.._ndim–1):                      *
 *        _x(i) = pbc(_x(i) + delta(i), side(i));                *
 *      dove side(i) è la lunghezza della scatola in quella       *
 *      dimensione, e pbc() applica le condizioni periodiche.    *
 *    - Aggiorna solo _x; non modifica _xold fino a che           *
 *      acceptmove() non viene chiamato.                         *
 *                                                              *
 * 4) void Particle::flip()                                      *
 *    - Inverte il valore dello spin: _spin = –_spin.            *
 *      Utilizzato nel modello di Ising (flip di un singolo spin).*
 *    - Non restituisce alcun valore.                            *
 *                                                              *
 * 5) void Particle::moveback()                                  *
 *    - Ripristina la posizione corrente alle coordinate “vecchie”*
 *      memorizzate in _xold: _x = _xold.                        *
 *    - Usato per annullare una proposta di spostamento MC LJ se  *
 *      la mossa viene rifiutata (metro() restituirà false).     *
 *                                                              *
 * 6) void Particle::acceptmove()                                *
 *    - Aggiorna _xold = _x, salvando la nuova posizione come    *
 *      “posizione precedente” per il prossimo passo.            *
 *    - Utile sia in MD (dopo Verlet, per fissare x_old = x_new) *
 *      sia in MC (dopo un move() accettato).                     *
 *                                                              *
 * 7) int Particle::getspin()                                    *
 *    - Restituisce il valore corrente di _spin (+1 o –1).       *
 *    - Utile per calcoli Ising.                                  *
 *                                                              *
 * 8) void Particle::setspin(int spin)                           *
 *    - Imposta esplicitamente lo spin a +1 o –1.                *
 *    - Void, non ritorna nulla.                                  *
 *                                                              *
 * 9) double Particle::getposition(int dim, bool xnew)           *
 *    - Se xnew == true, restituisce _x(dim), ovvero la          *
 *      coordinata corrente nella dimensione dim.                *
 *    - Se xnew == false, restituisce _xold(dim), cioè la        *
 *      posizione “precedente” nella dimensione dim.              *
 *    - Consente di accedere in lettura a entrambe le posizioni.  *
 *                                                              *
 * 10) void Particle::setposition(int dim, double position)      *
 *    - Imposta direttamente la coordinata corrente _x(dim)       *
 *      al valore “position”.                                    *
 *    - Non tocca _xold finché non viene chiamato acceptmove().    *
 *                                                              *
 * 11) void Particle::setpositold(int dim, double position)      *
 *    - Imposta la coordinata “vecchia” _xold(dim) su                 *
 *      “position”.                                               *
 *    - Utilizzato all’inizializzazione delle velocità in MD,    *
 *      per costruire x_old = x(t) – v·Δt.                         *
 *                                                              *
 * 12) double Particle::getvelocity(int dim)                     *
 *    - Restituisce la componente di velocità _v(dim).           *
 *    - Utile per calcoli di energia cinetica o istogrammi.       *
 *                                                              *
 * 13) vec Particle::getvelocity()                               *
 *    - Restituisce l’intero vettore di velocità _v (arma::vec    *
 *      di lunghezza _ndim).                                      *
 *    - Comodo per sommare ~½m v² in System::measure().           *
 *                                                              *
 * 14) void Particle::setvelocity(int dim, double velocity)      *
 *    - Imposta la componente di velocità _v(dim) al valore dato. *
 *    - Utilizzato sia all’inizializzazione delle velocità MD     *
 *      sia durante integrazione Verlet.                          *
 *                                                              *
 * 15) double Particle::pbc(double position, double side)        *
 *    - Applica le condizioni periodiche per una coordinata scalare:*
 *      Restituisce position – side * rint(position/side).         *
 *    - Riporta la particella dentro la scatola [0, side) o     *
 *      [−side/2, side/2] a seconda del sistema di riferimento.    *
 *    - Chiamato da translate() e da Verlet tramite System.       *
 *                                                              *
 * 16) double Particle::getvelocity_module()                     *
 *    - Calcola e restituisce il modulo della velocità vettoriale:*
 *        velocity_module = sqrt(Σ_{i=0.._ndim−1} v(i)²).          *
 *    - Utile per costruire l’istogramma di p(v) nel data blocking.*
 *                                                              *
 * 17) Conclusione                                              *
 *    - Tutti i metodi sono void o ritornano un dato scalare/int.*
 *    - La classe Particle espone operazioni di base su ciascuna *
 *      particella, quali spostamento, aggiornamento spin,        *
 *      gestione delle coordinate vecchie per Verlet/MC e accesso  *
 *      alle velocità.                                            *
 ****************************************************************/

#include <iostream>
#include <math.h>
#include "particle.h"

using namespace std;

void Particle :: initialize(){
   _spin = 1;
   _x.resize(_ndim);
   _xold.resize(_ndim);
   _v.resize(_ndim);
   return;
}

void Particle :: translate(vec delta, vec side){
   for(unsigned int i=0; i<_ndim; i++){
     _x(i) = pbc(_x(i) + delta(i), side(i));
   }
}

void Particle :: flip(){
   _spin = -1*this->getspin();
}

void Particle :: moveback(){
   _x = _xold;
}

void Particle :: acceptmove(){
   _xold = _x;
}

int Particle :: getspin(){
   return _spin;
}

void Particle :: setspin(int spin){
   _spin = spin;
   return;
}

double Particle :: getposition(int dim, bool xnew){
   if(xnew) return _x(dim);
   else return _xold(dim);
}

void Particle :: setposition(int dim, double position){
   _x(dim) = position;
   return;
}

void Particle :: setpositold(int dim, double position){
   _xold(dim) = position;
   return;
}

double Particle :: getvelocity(int dim){
   return _v(dim);
}

vec Particle :: getvelocity(){
   return _v;
}

void Particle :: setvelocity(int dim, double velocity){
   _v(dim) = velocity;
   return;
}

double Particle :: pbc(double position, double side){
  return position - side * rint(position / side);
}
double Particle :: getvelocity_module(){
   double velocity_module = 0;
   for(int i = 0; i < _ndim; i++) velocity_module += pow(_v(i),2);
   return sqrt(velocity_module);
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
