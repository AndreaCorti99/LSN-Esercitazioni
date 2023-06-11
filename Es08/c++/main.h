/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __QuantumMC__
#define __QuantumMC__

#include "random.h"
#include <vector>

using namespace std;


int seed[4];
Random rnd;

//medie a blocchi
int n_samples=1e5; //# campionamenti di energia
int n_blk=100; //blocchi di energia
int L=n_samples/n_blk;
vector<double> r(n_samples);//posizioni
vector<double> H(n_samples);
vector<double> H_sum_prog(n_blk);
vector<double> H_err_prog(n_blk);

//Simulated Annealing
const int n_step=1500; //numero di temperature per il SA
double H_SA[n_step];
double H_err_SA[n_step];
double Lmi=1.5, Lsigma=1.5;
double delta_mi, delta_sigma;
double mi[n_step], sigma[n_step];

//Parametri
double x=1.5; //inizializzazione casuale
int Acceptance=0;
double temp=1., beta=1./temp, Delta_beta=2.;

//pigreco
const double pi=3.1415927;

//functions
double Psi(double,double,double);
double DDPsi(double,double,double);
double Potential(double);
double H_Psi_over_Psi(double,double,double);
double H_Psi_over_PsiII(double,double,double);
void Rnd_Inizialization(void);
double error(vector<double>&,vector<double>&, int);
void MeanAndErr(vector<double>&,vector<double>&,vector<double>&,int,int);
bool MetrGauss(double&,double,Random&,double,double);
void Reset(void);
void Calculate_H(double,double, int, int);
int Minimum_Index(double[],int);

#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
