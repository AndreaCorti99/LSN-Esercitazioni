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
#include "mpi.h"


using namespace std;

void Rnd_Inizialization(Random& random, int index_line);

int dim =50;
int n_chromosome=3000;
int type=2;





int main (int argc, char *argv[]){

   int size, rank;

   //INIZIO PARALLELIZZAZIONE
   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD, &size);//ottengo num tot di processi
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);//ogni processo ottiene il proprio rank
   MPI_Status stat1, stat2;

   Random rnd;
   Rnd_Inizialization(rnd, rank);

   TSP ProblemSet(dim, rnd, type); //Inizializzo il problema, punti sulla circonferenza

   ProblemSet.LoadCities("American_capitals.dat");

   Population Pop(n_chromosome, dim, rnd); //Inizializzo una popolazione di cromosomi
cout << "chromosome 1 "; Pop.printChromosome_i(0); cout << "      rank " << rank << endl;
   Pop.Set_parameters(3.0, 0.5, 1*0.15, 1*0.1, 1*0.1, 1*0.1); //Imposto i parametri di evoluzione del sistema
   
   int Nmigr = 5; //definisce ogni quanto bisogna scambiare cromosomi tra popolazioni parallele

   vector<int>imesg (dim);
   vector<int>imesg2 (dim);
   int itag=1; int itag2=2;
   vector <int> which_swap = {0, 1, 2, 3};

   //EVOLUZIONE
   for (int i=1; i<=200; i++){
      //cout << "Generation: " << i << endl;
      Pop.Evolve(rnd,Pop,ProblemSet);

      if (i% Nmigr == 0) { //Ogni Nmigr generazioni si fa scambio tra continenti
      cout << "Generation: " << i << endl;	
         random_shuffle(which_swap.begin(), which_swap.end()); 

	 //scambio i migliori delle prime due popolazioni: invece che fare due scambi non si può fare Comm_split?
	 for(int j=0; j<dim; j++){	
	    imesg[j] = Pop.Get_i_Chrom(0).Get_i_gen(j);
	    imesg2[j] = Pop.Get_i_Chrom(0).Get_i_gen(j);
	 }
		
         if(rank==which_swap[1]){
	    MPI_Send(&imesg[0],dim,MPI_INTEGER,which_swap[0],itag,MPI_COMM_WORLD);
            MPI_Recv(&imesg2[0],dim,MPI_INTEGER,which_swap[0],itag2, MPI_COMM_WORLD,&stat2);
				//cout<<"messaggio1 = "<<imesg2[0]<<endl;
	 }
	 else if(rank==which_swap[0]){
	    MPI_Send(&imesg2[0],dim,MPI_INTEGER,which_swap[1],itag2, MPI_COMM_WORLD);
	    MPI_Recv(&imesg[0],dim,MPI_INTEGER,which_swap[1],itag, MPI_COMM_WORLD,&stat1);
				//cout<<"messaggio = "<<imesg[0]<<endl;
	 }
	
	 if(rank==which_swap[1]){ 
            Chromosome Chrom_back(dim);
            Chrom_back.SetGenes(imesg2);
            Pop.Set_i_Chrom(Chrom_back, 0); 
         }
	 else if(rank==which_swap[0]){ 
            Chromosome Chrom_back(dim);
            Chrom_back.SetGenes(imesg);
            Pop.Set_i_Chrom(Chrom_back, 0); 
         }

         Pop.MergeSort(Pop, ProblemSet, 0, n_chromosome - 1);
      }		

      PrintCities(ProblemSet, Pop, i, rank);
      PrintBestFitness(ProblemSet, Pop, i, rank);
      PrintBestAvFitness(ProblemSet, Pop, i, rank);
   }

 
   rnd.SaveSeed(); 

   MPI_Finalize(); //FINE PARALLELIZZAZIONE



   return 0;
}



//----------------------------------------------------------------//

//----------------------------------------------------------------//

//RANDOM NUMBER GENERATOR INIZIALIZATOR
void Rnd_Inizialization( Random & random, int index_line ) {

   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      for(int i=0; i<index_line; i++) {
         Primes >> p1 >> p2 ;
      }
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
