/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm> //serve per sort
#include "random.h"
#include "stat.h"

using namespace std;

//FUNZIONI
double S(double S0, double m, double s, double t, double W) {
   return S0*exp((m-s*s/2.)*t+s*W);
}

double Max(double a, double b) {
   if (a>b) {return a;}
   else return b;
}

void InitializeRand(Random& rnd); //Implementata sotto



 
int main (int argc, char *argv[]){

   Random rnd; 
   InitializeRand(rnd);


   //ESERCIZIO 3.1.1

   int M = 1e5;              //Total number of throws
   int NBlocks = 1e2;                 // Number of blocks
   int L = M/NBlocks;    		//# of numbers in a block

   //Parametri caratteristici
   double Mean = 0.1;
   double Sigma = 0.25;
   double S0 = 100.;
   double T = 1.;
   double k = 100.;

   Stat st;

   vector<double> W(M);
   for (int i=0; i<M; i++) {
      W[i] = rnd.Gauss(0,T);
   }   

   vector<double> Spesa(M);
   for (int i=0; i<M; i++) {
      Spesa[i] = S(S0,Mean,Sigma,T,W[i]);
   }   

   vector<double> C(M);
   vector<double> P(M);
   for (int i=0; i<M; i++) {
      C[i] = exp(-Mean*T)*Max(0.,Spesa[i]-k);
      P[i] = exp(-Mean*T)*Max(0.,k-Spesa[i]);
   }  

   vector<double> MeanC(NBlocks);
   vector<double> MeanP(NBlocks);
   vector<double> ErrorC(NBlocks);
   vector<double> ErrorP(NBlocks);

   st.MeanAndErr(C,MeanC,ErrorC,NBlocks,L);
   st.MeanAndErr(P,MeanP,ErrorP,NBlocks,L);

   ofstream outfile311("outfile311.txt");
    for (int i = 0; i < NBlocks; ++i) {
       outfile311 << MeanC[i] << "\t" << ErrorC[i] << "\t" << MeanP[i] << "\t" << ErrorP[i] << endl;
    }
    outfile311.close();


   //ESERCIZIO 3.1.2

   int NIntervalli = 100;
   double Incremento = T/NIntervalli;
   
   vector<double> SpesaD(M);
   for (int i=0; i<M; i++) {
      SpesaD[i] = S0;
   }      

   for (int i=0; i<M; i++) {
      double Zi = rnd.Gauss(0,1);
      for (int j=0; j<NIntervalli; j++) {
         SpesaD[i] = SpesaD[i]*exp((Mean-1./2.*Sigma*Sigma)*Incremento + Sigma*Zi*Incremento);
      }   
   }   

   for (int i=0; i<M; i++) {
      C[i] = exp(-Mean*T)*Max(0.,SpesaD[i]-k);
      P[i] = exp(-Mean*T)*Max(0.,k-SpesaD[i]);
   }  

   fill(MeanC.begin(), MeanC.end(), 0.0);
   fill(MeanP.begin(), MeanP.end(), 0.0);
   fill(ErrorC.begin(), ErrorC.end(), 0.0);
   fill(ErrorP.begin(), ErrorP.end(), 0.0);

   st.MeanAndErr(C,MeanC,ErrorC,NBlocks,L);
   st.MeanAndErr(P,MeanP,ErrorP,NBlocks,L);

   ofstream outfile312("outfile312.txt");
    for (int i = 0; i < NBlocks; ++i) {
       outfile312 << MeanC[i] << "\t" << ErrorC[i] << "\t" << MeanP[i] << "\t" << ErrorP[i] << endl;
    }
    outfile312.close();



   rnd.SaveSeed(); /*cambia i semi a fine lavoro*/

   return 0;
}




// ------------------------------//
// --------- FUNZIONI -----------//
// ------------------------------//

void InitializeRand(Random & rnd) {

   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ; //li legge dal file
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property; //variabile stringa chiamata property
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];//dai i semi all'algoritmo, presi dal file "seed.in"
            rnd.SetRandom(seed,p1,p2); 
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
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
