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
 * random.cpp: Implementazione di un generatore di numeri         *
 * casuali basato su MLCG con metodo di trasmissione del trasporto *
 *                                                              *
 * Questo file definisce la classe Random, che fornisce          *
 * funzionalità per generare numeri casuali uniformi e           *
 * gaussiani, salvare lo stato del generatore e inizializzarlo   *
 * con semi esterni.                                              *
 *                                                              *
 * Struttura del file:                                           *
 *                                                              *
 * 1) Inclusioni                                                *
 *    - <iostream>, <fstream>, <cmath>, <cstdlib>: librerie       *
 *      standard per I/O, funzioni matematiche e EXIT_FAILURE.    *
 *    - "random.h": dichiarazione della classe Random.            *
 *                                                              *
 * 2) Costruttore e Distruttore                                  *
 *    - Random(): costruttore vuoto che non fa nulla di specifico.*
 *    - ~Random(): distruttore vuoto.                             *
 *                                                              *
 * 3) void Random::SaveSeed()                                    *
 *    - Apre il file ../OUTPUT/seed.out in modalità scrittura.    *
 *    - Se aperto correttamente, scrive i quattro parametri        *
 *      l1, l2, l3, l4, che rappresentano lo stato interno del RNG.*
 *    - Altrimenti, stampa un messaggio d’errore.                 *
 *    - Chiude il file.                                           *
 *                                                              *
 * 4) double Random::Gauss(double mean, double sigma)            *
 *    - Genera due numeri uniformi s, t ∈ (0,1) tramite Rannyu(). *
 *    - Applica la trasformazione Box-Muller:                     *
 *        x = sqrt(-2 log(1–s)) cos(2π t).                        *
 *    - Restituisce mean + x·sigma, generando così una distribuzione*
 *      gaussiana N(mean, sigma²).                                 *
 *                                                              *
 * 5) double Random::Rannyu(double min, double max)              *
 *    - Restituisce un numero uniforme in [min, max] calcolando:  *
 *        min + (max–min) * Rannyu().                              *
 *                                                              *
 * 6) double Random::Rannyu(void)                                 *
 *    - Implementazione del generatore lineare congruenziale       *
 *      a quattro parametri con passo di trasporto:              *
 *        • Calcola i1,i2,i3,i4 come combinazioni lineari di       *
 *          semi l1..l4 e le componenti di trasporto n1..n4.       *
 *        • Applica modulo 4096 su ciascun aggiornamento di l1..l4.*
 *        • Costruisce r = 2^(-12) [l1 + 2^(-12)(l2 + 2^(-12)(l3 + 2^(-12) l4))].*
 *        • Restituisce r ∈ [0,1).                                 *
 *                                                              *
 * 7) void Random::SetRandom(int *s, int p1, int p2)              *
 *    - Imposta i parametri del generatore:                        *
 *        • m1=502, m2=1521, m3=4071, m4=2107 (moduli interni).    *
 *        • l1..l4 = s[0..3], dove s è un array di quattro semi.   *
 *        • n1=0, n2=0, n3=p1, n4=p2 (semi di trasporto).           *
 *    - Non restituisce nulla (void).                              *
 *                                                              *
 * 8) Macro di chiusura                                            *
 *    - fine di #ifdef __Random__ e #endif                         *
 ****************************************************************/
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "random.h"

using namespace std;

Random :: Random(){}

Random :: ~Random(){}

void Random :: SaveSeed(){
   ofstream WriteSeed;
   WriteSeed.open("../OUTPUT/seed.out");
   if (WriteSeed.is_open()){
      WriteSeed << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
   } else cerr << "PROBLEM: Unable to open random.out" << endl;
  WriteSeed.close();
  return;
}

double Random :: Gauss(double mean, double sigma) {
   double s=Rannyu();
   double t=Rannyu();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   return mean + x * sigma;
}

double Random :: Rannyu(double min, double max){
   return min+(max-min)*Rannyu();
}

double Random :: Rannyu(void){
  const double twom12=0.000244140625;
  int i1,i2,i3,i4;
  double r;

  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;
  l4 = i4%4096;
  i3 = i3 + i4/4096;
  l3 = i3%4096;
  i2 = i2 + i3/4096;
  l2 = i2%4096;
  l1 = (i1 + i2/4096)%4096;
  r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));

  return r;
}

void Random :: SetRandom(int * s, int p1, int p2){
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0];
  l2 = s[1];
  l3 = s[2];
  l4 = s[3];
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
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
