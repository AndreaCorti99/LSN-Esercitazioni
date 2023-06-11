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
#include <cmath>
#include <functional> //serve per definire funzione di funzione
#include <iomanip> //serve per setprecision
#include "TSP.h"


using namespace std;

void Rnd_Inizialization(Random& random, int* seed);


int main (/*int argc, char *argv[]*/){

   //Inizializzazione generatore di numeri pseudocausali
   int seed_[4];
   Random rnd;
   Rnd_Inizialization(rnd, seed_);

   //Definizione dei parametri e del tipo del problema tramite terminale
   int dim = 34;

   int n_chromosome;
   cout << "Insert the desired Number of Chromosomes: ";
   cin >> n_chromosome;
   //int n_chromosome=300;
   string Chr_dim = to_string(n_chromosome);

   int type;
   cout << "Select the problem type: 0 <------> Euclidean Circle " << endl;
   cout << "                         2 <------> Euclidean Square " << endl;
   cin >> type;
   //int problem_type=0;
   string Problem_type; 
   if (type==0) {Problem_type = "Circle";}
   if (type==2) {Problem_type = "Square";}


   TSP ProblemSet(dim, rnd, type); //Inizializzo il problema, punti sulla circonferenza

   Population Pop(n_chromosome, dim, rnd); //Inizializzo una popolazione di cromosomi
   Pop.Set_parameters(2.0, 0.5, 0.15, 1*0.1, 1*0.1, 1*0.1); //Imposto i parametri di evoluzione del sistema
   

   //EVOLUZIONE
   for (int i=1; i<=300; i++){
      cout << "Generation: " << i << endl;
      Pop.Evolve(rnd,Pop,ProblemSet);
      PrintCities(ProblemSet, Pop, i, Problem_type, Chr_dim);
      PrintBestFitness(ProblemSet, Pop, i, Problem_type, Chr_dim);
      PrintBestAvFitness(ProblemSet, Pop, i, Problem_type, Chr_dim);
   }


   rnd.SaveSeed(); /*cambia i semi a fine lavoro*/

   return 0;
}



//----------------------------------------------------------------//

//----------------------------------------------------------------//

//RANDOM NUMBER GENERATOR INIZIALIZATOR
void Rnd_Inizialization( Random & random, int* seed ) {

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
            random.SetRandom(seed,p1,p2); 
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
