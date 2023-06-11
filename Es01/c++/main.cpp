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
#include <iomanip>
#include <algorithm> //serve per sort
#include "random.h"
#include "stat.h"

using namespace std;

void InitializeRand(Random& rnd); //Implementata sotto

 
int main (int argc, char *argv[]){

   Random rnd; /*dichiari la variabile rnd appartenente alla classe random, presa dai file random.cpp e random.h */
   InitializeRand(rnd);


   //INIZIO dell'esercizio vero e proprio


   //ESERCIZIO 1.1.1

   int M = 1e6;              //Total number of throws
   int N = 1e2;                 // Number of blocks
   int L = M/N;    		//# of numbers in a block

   //contiene i dati casuali grezzi
   vector<double> r(M);
   for (int i = 0; i < M; i++) {
        r[i] = rnd.Rannyu(); // U[0,1) uniform distribution
   }

   //tengo salvati nel main solo i vettori contenenti i dati che poi verranno salvati su file. Gli altri creati e distrutti all'interno della funzione MeanAndErr	
   vector<double> sum_prog(N);
   vector<double> err_prog(N);
   Stat st;
   st.MeanAndErr(r,sum_prog,err_prog,N,L);

    ofstream outfile1("output/outfileMeanValue.txt");
    for (int i = 0; i < N; ++i) {
        outfile1 << sum_prog[i] << setw(20) << err_prog[i] << endl;
    }
    outfile1.close();


   //ESERCIZIO 1.1.2

   fill(sum_prog.begin(), sum_prog.end(), 0.0);
   fill(err_prog.begin(), err_prog.end(), 0.0);

   vector<double> s(M);
   for(int i=0; i<M; i++) {
      s[i] = pow(r[i]-0.5,2);
   }

   st.MeanAndErr(s, sum_prog, err_prog, N, L);

   ofstream outfileSigma1("output/outfileSigma.txt");
   for (int i = 0; i < N; ++i) {
       outfileSigma1 << sum_prog[i] << setw(20) << err_prog[i] << endl;
   }
   outfileSigma1.close();


   //ESERCIZIO 1.1.3

   int Nchi=1000; //numero di chi quadrati che si vogliono calcolare, ricorda che se si vuole statistica sensata i dati da usare per ogni chi devono essere molti di più degli intervalli disponibili. Per esempio se abbiamo 1e6 dati generati, e facciamo 1e3 chi quadri, ci rimangono 1e3 dati per ogni set e ci aspettiamo 10 dati circa a intervallo.
   int Lnew = M/Nchi;	
   vector<vector<double>> matrix(N, vector<double>(Nchi, 0.0));//crea una matrice NxNchi (come vettore di vettori) piena di zeri

   for (int j = 0; j < Nchi; j++) {
      vector<double> VtoSort(r.begin() + j * Lnew, r.begin() + (j + 1) * Lnew);
      sort(VtoSort.begin(), VtoSort.end());
      for (int i = 0; i < N; i++) {
         double count = 0;
         for (int k = 0; k < Lnew; k++) {
            if ((i / double(N) <= VtoSort[k]) && (VtoSort[k] < (i + 1) / double(N))) {
                count++;
            }
         }
         matrix[i][j] = count;
      }
   }

   vector<double> ChiVect(Nchi);
   for (int j = 0; j < Nchi; j++) {
      double sum = 0.;
      for (int i = 0; i < N; i++) {
          sum += pow(matrix[i][j] - double(Lnew/N), 2) / double(double(Lnew)/double(N));
      }
      ChiVect[j] = sum;
   }

   ofstream outfileChi2("output/outfileChi2.txt");
   for (int i = 0; i < Nchi; ++i) {
        outfileChi2 << ChiVect[i] << endl;
   }
   outfileChi2.close();



   //ESERCIZIO 1.2.1

   int Nrep = 1e5;
   int dimension[] = {1, 2, 10, 100};
   double Mean = 1.;

   for(int l=0; l<4; l++) {

      vector<double> S_N(Nrep);
      for (int i=0; i<Nrep; i++) {
         for (int j=0; j<dimension[l]; j++) {
            S_N[i]+=rnd.Exp(Mean);
         }
         S_N[i]=S_N[i]/dimension[l];
      }

      string Index = to_string(dimension[l]);
      ofstream outfileLCTexp("output/ExpLCT/outfileLCTexp"+Index+".txt");
      for (int i = 0; i < Nrep; ++i) {
           outfileLCTexp << S_N[i] << endl;
      }
      outfileLCTexp.close();
   }


   double G = 1.;
   for(int l=0; l<4; l++) {

      vector<double> S_N(Nrep);
      for (int i=0; i<Nrep; i++) {
         for (int j=0; j<dimension[l]; j++) {
            S_N[i]+=rnd.Lor(G);
         }
         S_N[i]=S_N[i]/dimension[l];
      }

      string Index = to_string(dimension[l]);
      ofstream outfileLCTlor("output/LorLCT/outfileLCTlor"+Index+".txt");
      for (int i = 0; i < Nrep; ++i) {
           outfileLCTlor << S_N[i] << endl;
      }
      outfileLCTlor.close();
   }


   //ESERCIZIO 1.2.2 

   double Length = 1.;
   double d = 2.;
   M = 1e5;         //Total number of throws
   N = 1e2;         // Number of blocks
   L = M/N; 

   vector<double> b(M);
   for (int i = 0; i < M; i++) {
        b[i] = 1.*d*rnd.Rannyu(); // U[0,1) uniform distribution
   }

   vector<double> l(M);

   for (int i = 0; i < M; i++) {
      l[i] = Length*sin(2.*rnd.UnPhiAR());
   }

   vector<double> ave_(N);
   vector<double> av2_(N);
   vector<double> sum_prog_(N); //nuovi vettori
   vector<double> su2_prog_(N);
   vector<double> err_prog_(N);

    for (int i = 0; i < N; i++) {
        int Nhit = 0;
        for (int j = 0; j < L; j++) {
            int k = j + i * L;
            if ( b[k]+l[k] < 0. || b[k]+l[k] > d ) {
                Nhit++;
            }
        }
        ave_[i] = 2.0 * Length * L / (d * Nhit); // calculation of pi for the i-th block
        av2_[i] = pow(ave_[i], 2); // (r_i)^2
    }

   for (int i = 0; i < N; i++) {
        for (int j = 0; j <= i; j++) {
            sum_prog_[i] += ave_[j];
            su2_prog_[i] += av2_[j];
        }
        sum_prog_[i] /=(i+1);
        su2_prog_[i] /=(i+1);
        err_prog_[i] = st.error(sum_prog_, su2_prog_, i);
    }

    ofstream outfile13("output/outfile13.txt");
    for (int i = 0; i < N; ++i) {
        outfile13 << sum_prog_[i] << "\t" << err_prog_[i] << endl;

    }
    outfile13.close();


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
