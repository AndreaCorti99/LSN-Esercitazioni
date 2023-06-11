#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "stat.h"

using namespace std;

// Default constructor, does not perform any action
Stat :: Stat(){}
// Default destructor, does not perform any action
Stat :: ~Stat(){}

//Calcola l'errore su media a blocchi
double Stat :: error(vector<double> &AV, vector<double> &AV2, int n) {
   if (n == 0) {
        return 0.0;
   } else {
        return sqrt((AV2[n] - AV[n]*AV[n])/n); // il /n dovrebbe corrispondere al /N-1
   }
}

//Aggiorna gli array passati come argomenti con medie e errori, partendo dal vettore di numeri casuali r
void Stat:: MeanAndErr(vector<double> &r, vector<double> &Mean, vector<double> &Errors, int N, int L) {

   Stat st;
   vector<double> ave(N);
   vector<double> av2(N);
   vector<double> su2_prog(N);

   for (int i = 0; i < N; i++) {
      double sum1 = 0.0;
      for (int j = 0; j < L; j++) {
         int k = j + i * L;
         sum1 += r[k];
      }
      ave[i] = sum1 / L;
      av2[i] = pow(ave[i], 2);
   }

   for (int i = 0; i < N; i++) {
        for (int j = 0; j <= i; j++) {
            Mean[i] += ave[j];
            su2_prog[i] += av2[j];
        }
        Mean[i]/=(i+1); //realizza il /n nel calcolo della media sui valori dei blocchi
        su2_prog[i]/=(i+1);
        Errors[i] = st.error(Mean, su2_prog, i);
    }
    return;
}

   
