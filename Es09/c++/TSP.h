/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __TSP__
#define __TSP__

#include <iostream>
#include "random.h"
#include <vector>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <algorithm>

using namespace std;

//funzione Periodic Boundary Conditions
int Pbc(int pos, int Period){
   if(pos>=0)    pos = pos%(Period); 
   if(pos < 0 && pos%(Period)!=0)      pos = Period + pos%Period;
   if(pos < 0 && pos%(Period)==0)      pos = pos%Period;
   return pos;
}



//----------------------------------------//
//---------- Classe pointCoord -----------//
//----------------------------------------//

//Rappresenta un generico punto in due dimensioni

class pointCoord {

public:

   pointCoord(double x, double y): m_x(x), m_y(y) {};
   pointCoord(double theta) : m_theta(theta) {}; //può essere comodo quando si considerano punti su una circonferenza di raggio fissato
   pointCoord(): m_x(0), m_y(0) {};
   pointCoord(const pointCoord& p): m_x(p.GetX()), m_y(p.GetY()) {};
   ~pointCoord(){};
		
   double distance ( pointCoord p) const {return sqrt( pow(m_x-p.GetX(),2) + pow(m_y-p.GetY(),2) );}; //calcola la distanza euclidea tra due punti nel piano
   double distanceRad (pointCoord p) const {return m_theta-p.GetTheta();}; //calcola la distanza in radianti tra due punti sulla circonferenza
   double Cart_to_Rad (double x, double y) const {return atan2(y, x);};
		
   double GetX() const {return m_x;};
   double GetY() const {return m_y;};
   double GetTheta() const {return m_theta;};
		
   void SetX(double x) {m_x=x;};
   void SetY(double y) {m_y=y;};
   void SetTheta(double theta) {m_theta=theta;};
		
protected:
	
   double m_x, m_y;
   double m_theta;

};



//----------------------------------------//
//----------- Classe Chromosome ----------//
//----------------------------------------//

//Classe Cromosoma, contiene metodi generali che caratterizzano qualsiasi algoritmo genetico.

class Chromosome {

public:

    Chromosome() {}; //Necessario perché quando inizializzo Population, prima di chiamare il costruttore sotto, chiama questo di default

    Chromosome(int N_gen) { //Costruttore
        n_genes = N_gen;
        genes.resize(N_gen); //genes is a double<vector>. Resize establish its dimension dinamically
        initializeGenes();
    }

    ~Chromosome() {} //Distruttore

    void initializeGenes() { //Inizializza il vettore genes
        for (int i = 0; i < n_genes; i++) {
            genes[i] = i; // Initialize genes with random integer values between 0 and N_gen
        }
    }

    //Accesso e Assegnazione per le variabili genes e n_genes
    vector<int> GetGenes() const { return genes; }
    int Get_i_gen(int i)  const {return genes[i];}
    int GetN_gen() const { return n_genes; }
    void SetGenes( vector<int> gene_to_set ) { genes = gene_to_set;}

    //CHECK
    void Check() const {
       if(genes[0]!=0) { cout << "First element moved!" << endl;}
       for(int i=1; i<n_genes; i++) {
          if( count(genes.begin()+1, genes.end(), i != 1 )) { cout << "Repetition or lack of elements!";}
       }
    }

    //MUTAZIONI E SCAMBI

    //1. Permuta due geni con posizioni a caso
    void Permutation (Random & rnd) {

       int index1 = floor(rnd.Rannyu(1.,n_genes)); //la prima città non deve essere mossa, allora parto da 1. e non da 0.
       int index2 = floor(rnd.Rannyu(1.,n_genes));
       int Back_variable = genes[index1];
       genes[index1]=genes[index2];
       genes[index2]=Back_variable;
    }

    //1.bis Permuta due geni con una posizione fissata
    void Permutation (Random & rnd, int m) {

       int n = floor(rnd.Rannyu(1.,n_genes)); //la prima città non deve essere mossa, allora parto da 1. e non da 0.
       int Back_variable = genes[m];
       genes[m]=genes[n];
       genes[n]=Back_variable;
    }

    //2. Trasla di +n posizioni m geni contigui
    void Shift (Random & rnd, int n, int m) {

       //int n = floor(rnd.Rannyu(1.,n_genes-1)); //scommenta se vuoi la funzione senza int inputs
       //int m = floor(rnd.Rannyu(1.,n_genes-1)); //interi tra 1, 2, ..., n_genes-2

       int i_start = floor(rnd.Rannyu(1.,n_genes)); //indice da cui far partire lo shift, va da 1, 2, ..., n_genes-1
  
       vector<int> genes_cut(n_genes-1);
       vector<int> temp(n_genes - 1); // Vettore temporaneo di dimensione minore
       i_start=i_start-1; //nel vettore tagliato l'indice diminuisce di una unità

       for(int i = 1; i<n_genes; i++) { //riempio nuovo vettore con la posizione iniziale tagliata
          genes_cut[i-1] = genes[i];
       }

       for (int j=0; j<n; j++) {//calcolo ricorsivamente ad ogni singolo shift la posizione in genes

          for(int i = 0; i<m; i++) { //sposto a destra il blocco di una posizione
             temp[Pbc(i_start+i+1, n_genes-1)] = genes_cut[Pbc(i_start+i , n_genes-1)];
          }

          temp[Pbc(i_start, n_genes-1)] = genes_cut[Pbc(i_start+m , n_genes-1)]; //sposto prima del blocco il termine che stava dopo

          for(int i = 0; i<n_genes-2-m; i++) { //lascio fermi gli altri termini
             temp[Pbc(i_start+m+1+i, n_genes-1)] = genes_cut[Pbc(i_start+m+1+i, n_genes-1)]; 
          }

          for(int i = 1; i<n_genes; i++) { //riempio nuovo vettore con la posizione iniziale tagliata
             genes_cut[i-1] = temp[i-1];
          }

          i_start++;
       }

       for (int i = 1; i < n_genes; i++) { //Rimetti il vettore di supporto nel vettore originale
          genes[i] = temp[i - 1];
       }
    } 

    //3. Permuta m geni contigui
    void Permutation_m (Random & rnd) {
       int m = floor(rnd.Rannyu(1.,n_genes-1)); 
       Shift(rnd,m,m);
    }

    //4. Inverti un blocco di m geni
    void Inversion (Random & rnd) {

       int m = floor(rnd.Rannyu(1.,n_genes-1)); 
       int i_start = floor(rnd.Rannyu(1.,n_genes)); //indice da cui far partire lo shift

       vector<int> genes_cut(n_genes-1);
       vector<int> temp(n_genes - 1); // Vettore temporaneo di dimensione minore
       i_start=i_start-1; //nel vettore tagliato l'indice diminuisce di una unità

       for(int i = 1; i<n_genes; i++) { //riempio nuovo vettore con la posizione iniziale tagliata
          genes_cut[i-1] = genes[i];
       }

       for (int i = i_start; i < i_start+m; i++) { //ciclo sul blocco da invertire
             temp[Pbc(i_start+(m-1)-(i-i_start), n_genes-1)] = genes_cut[Pbc(i , n_genes-1)]; 
       }

       for (int i = i_start+m; i < i_start+n_genes-1; i++) { //lascio fermi quelli che non devono muoversi
             temp[Pbc(i,n_genes-1)] = genes_cut[Pbc(i,n_genes-1)]; 
       }

       for (int i = 1; i < n_genes; i++) { //Rimetti il vettore di supporto nel vettore originale
          genes[i] = temp[i - 1];
       }       
    }

    //5. Crossover
    void Crossover(Random & rnd, Chromosome & parent2) { //modificato sia il chr su cui è applicato il crossover, sia il chr argomento della funzione.

       int i_start = floor(rnd.Rannyu(1.,n_genes-1)); //indice da cui far partire il crossover
       int support_index=0;

       vector<int> parent1_genes = genes;//il primo "for" modifica genes, quindi creo copia che serve nel secondo 
       vector<int> parent2_genes = parent2.GetGenes();

       for(int j=1; j<n_genes; j++) { 
          if( count(genes.begin()+1, genes.begin()+i_start, parent2.Get_i_gen(j)) == 0 ) { //controllo che l'elemento j di parent non sia nella parte fissata di genes
             genes[i_start+support_index]=parent2.Get_i_gen(j); //in tal caso procedo a mettere secondo l'ordine in cui compaiono in parent
             support_index++;                
          }              
       }

       support_index=0;
       for(int j=1; j<n_genes; j++) {
          if( count(parent2_genes.begin()+1, parent2_genes.begin()+i_start, parent1_genes[j]) == 0 ) { 
             parent2_genes[i_start+support_index]=parent1_genes[j];
             support_index++;                
          }              
       }
       parent2.SetGenes(parent2_genes);
       
    }


private:

    int n_genes;
    vector<int> genes;

};



//----------------------------------------//
//------------- Classe TSP ---------------//
//----------------------------------------//

//Classe commesso viaggiatore, contiene metodi specifici caratteristici del problema considerato

class TSP {

public:

   TSP( int n_genes , Random & rnd, int type) { //costruttore
      cities.resize( n_genes );
      SetProblemType(type);
      initializeCities(rnd);
   };

   ~TSP() {}; //distruttore

   //Accesso alle variabili
   int GetProblemType () const {return problem_type;}
   vector<pointCoord> GetCities () const { return cities; }
   pointCoord GetCity (int index) const { return cities[index]; } 

   //Altri metodi 
   void SetProblemType( int type ) { //imposta il tipo di problema

      if (type!=0 && type!=1 && type!=2) {
         cout << "Select a valid TSP: 0 <------> Euclidean Circle " << endl;
         cout << "                    1 <------> Circle Manifold "  << endl;
         cout << "                    2 <------> Euclidean Square " << endl;
         return;
      }
      
      problem_type = type;
   }


   void initializeCities(Random & rnd) { //inizializza il vettore "cities" in base al problema scelto
      
      if(problem_type==0) {
         for (int i = 0; i < static_cast<int>(cities.size()); i++) {
             double phi = 2.*M_PI*rnd.Rannyu();
             cities[i].SetX(cos(phi));  //Assign cities with a random couple of coordinate on a unit circle. The initialization is done by default constructor
             cities[i].SetY(sin(phi));
         }
      }

      if(problem_type==1) {
         for (int i = 0; i < static_cast<int>(cities.size()); i++) {
             cities[i].SetTheta(2.*M_PI*rnd.Rannyu());  //Assign cities with a random couple of coordinate on a unit circle. Initialization done by default constructor
         }
      }

      if(problem_type==2) {
         for (int i = 0; i < static_cast<int>(cities.size()); i++) {
             cities[i].SetX(rnd.Rannyu());  //Assign cities with a random couple of coordinate on a unit circle. Initialization done by default constructor
             cities[i].SetY(rnd.Rannyu());
         }
      }

   }


   //Funzioni fitness
   double Fitness1(Chromosome Chrom) const { //funzione fitness L1

      double distance=0.;
      for(int i=0; i<Chrom.GetN_gen()-1; i++) {
         int m =  Chrom.Get_i_gen(i) ; //i-esima elemento del cromosoma, gli sarà associato un intero che specifica l'ordine in cui viene visitata la città
         int n =  Chrom.Get_i_gen(i+1) ;
         distance += cities[m].distance(cities[n]);
      }
      //ultimo step tiene conto delle condizioni al contorno
      int m =  Chrom.Get_i_gen(Chrom.GetN_gen()-1) ; 
      int n =  Chrom.Get_i_gen(0) ;
      distance += cities[m].distance(cities[n]);

      return distance;

   }; 

   double Fitness2(Chromosome Chrom) const { //funzione fitness L2

      double distance2=0.;
      for(int i=0; i<Chrom.GetN_gen()-1; i++) {
         int m =  Chrom.GetGenes()[i] ; //i-esima elemento del cromosoma, gli sarà associato un intero che specifica l'ordine in cui viene visitata la città
         int n =  Chrom.GetGenes()[i+1] ;
         distance2 += pow(cities[m].distance(cities[n]),2.);
      }
      //ultimo step tiene conto delle condizioni al contorno
      int m =  Chrom.GetGenes()[Chrom.GetN_gen()] ; 
      int n =  Chrom.GetGenes()[0] ;
      distance2 += pow(cities[m].distance(cities[n]),2.);

      return distance2;

   }; 

protected:

private:

  int problem_type;
  vector<pointCoord> cities;

};



//----------------------------------------//
//---------- Classe Popolazione ----------//
//----------------------------------------//

//Classe Popolazione, contiene metodi utili alla nascita e all'evoluzione di una popolazione di cromosomi

class Population {

public:

    Population(int n_chr, int N_gen, Random & rnd) { //Costruttore
       population.resize(n_chr); // Resize population vector

       for (int i = 0; i < n_chr; i++) {
          population[i] = Chromosome(N_gen); //Crea istanze ai cromosomi separate per ogni elemento 
       }

       for (int i = 0; i<n_chr; i++) { //disordina tutti i cromosomi, il loro costruttore li crea ordinati
          for(int j = 0; j<N_gen; j++) {
             population[i].Permutation(rnd);
          }
       }

       n_chromosomes = n_chr;
       n_genes = N_gen;
    }
    
    Population(const Population& other) {
        population = other.population; 
        n_chromosomes = other.n_chromosomes;
        n_genes = other.n_genes;
    }
    
    Population() {}
    ~Population() {} //Distruttore


    //ACCESSO ai membri della classe
    vector<Chromosome>& GetChromosomes() { return population; } //NON mettere const, la mutazione deve modificare il cromosoma e l'accesso avviene con questo metodo
    void SetPopulation(Population Pop_to_set) { population = Pop_to_set.population; }
    Chromosome Get_i_Chrom(int i) const { return population[i]; }
    void Set_i_Chrom(Chromosome Chr_to_set, int i) { population[i] = Chr_to_set; }
    void Set_parameters(double a, double b, double c, double d, double e, double f) { p=a; Pc=b; Pp=c; Ps=d; Ppm=e; Pi=f; }


    //ORDINA il vettore secondo fitness
    //NB Unico punto in cui la classe Population dipende da un metodo specifico del set del problema!!
    void Sort( TSP ProblemSet ) {
       for(int i=0; i<n_chromosomes; i++) {//metto i cromosomi più "in forma" ai primi indici della popolazione
          for(int j=i; j<n_chromosomes; j++) {
             if(ProblemSet.Fitness1(population[j]) < ProblemSet.Fitness1(population[i])) {//Scambio i membri della popolazione
                Chromosome Back = population[i];
                population[i] = population[j];
                population[j] = Back;
             } 
          }
       }  
    }  

    void Merge(Population& Pop, TSP& ProblemSet, int left, int mid, int right) {
       int i, j, k;
       int n1 = mid - left + 1;
       int n2 = right - mid;

       // Create temporary arrays
       Chromosome* leftArr = new Chromosome[n1];
       Chromosome* rightArr = new Chromosome[n2];

       // Copy data to temporary arrays
       for (i = 0; i < n1; i++)
           leftArr[i] = Pop.population[left + i];
       for (j = 0; j < n2; j++)
           rightArr[j] = Pop.population[mid + 1 + j];

       // Merge the temporary arrays back into Pop.population[left..right]
       i = 0;
       j = 0;
       k = left;
       while (i < n1 && j < n2) {
           // Compare Fitness1 values of the chromosomes
           if (ProblemSet.Fitness1(leftArr[i]) <= ProblemSet.Fitness1(rightArr[j])) {
               Pop.population[k] = leftArr[i];
               i++;
           } else {
               Pop.population[k] = rightArr[j];
               j++;
           }
           k++;
       }

       // Copy the remaining elements of leftArr[], if any
       while (i < n1) {
           Pop.population[k] = leftArr[i];
           i++;
           k++;
       }

       // Copy the remaining elements of rightArr[], if any
       while (j < n2) {
           Pop.population[k] = rightArr[j];
           j++;
           k++;
       }

       // Free the temporary arrays
       delete[] leftArr;
       delete[] rightArr;
   }

   void MergeSort(Population& Pop, TSP& ProblemSet, int left, int right) {
       if (left < right) {
           int mid = left + (right - left) / 2;

           // Sort first and second halves
           MergeSort(Pop, ProblemSet, left, mid);
           MergeSort(Pop, ProblemSet, mid + 1, right);

           // Merge the sorted halves
           Merge(Pop, ProblemSet, left, mid, right);
       }
   }



    //EVOLUZIONE della popolazione                            
    void Evolve( Random & rnd, Population & Pop, TSP ProblemSet) {

       Population offspring(Pop); //inizializzazione utile solo ad allocare memoria, poi le informazioi vengono sovrascritte
 
       int n_couple = int(n_chromosomes/2); 
       if(n_chromosomes%2 != 0) { n_couple= (n_chromosomes-1)/2; } 

       for(int i=0; i<n_couple; i++) { //ciclo di CROSSOVER, porta a n_chromosomes figli
      
          Population backup(Pop); //copia backup di Pop

          //SELEZIONE
          int m = int(n_chromosomes * pow(rnd.Rannyu(),p)) ; 
          int n = int(n_chromosomes * pow(rnd.Rannyu(),p)) ;

          if( rnd.Rannyu() < Pc ) { backup.Crossover(rnd,m,n); } //Probabilità di Crossover > 50%

          offspring.Set_i_Chrom( backup.Get_i_Chrom(m) , 2*i); //metto cromosomi generati (m,n) riordinati nella pop dei figli (posizioni crescenti,2i e 2i+1)
          offspring.Set_i_Chrom( backup.Get_i_Chrom(n) , 2*i+1);



          //MUTAZIONI sui figli
          //1. Permuta
          int j = floor(rnd.Rannyu(1.,n_genes));
          if( rnd.Rannyu() < Pp ) { offspring.population[2*i].Permutation(rnd, j); } //Permuto il gene con un altro a caso
          if( rnd.Rannyu() < Pp ) { offspring.population[2*i+1].Permutation(rnd, j); }

          //2. Trasla
          int s = floor(rnd.Rannyu(1.,n_genes-1)); //interi tra 1, 2, ..., n_genes-2
          int t = floor(rnd.Rannyu(1.,n_genes-1)); //
          if( rnd.Rannyu() < Ps ) { offspring.population[2*i].Shift(rnd, s, t); } //shifto di n posizioni un gruppo di m geni
          if( rnd.Rannyu() < Ps ) { offspring.population[2*i+1].Shift(rnd, s, t); }

          //3. Permuta un gruppo di geni
          if( rnd.Rannyu() < Ppm ) { offspring.population[2*i].Permutation_m(rnd); } 
          if( rnd.Rannyu() < Ppm ) { offspring.population[2*i+1].Permutation_m(rnd); }

          //4. Inverti un gruppo di geni
          if( rnd.Rannyu() < Pi ) { offspring.population[2*i].Inversion(rnd); }
          if( rnd.Rannyu() < Pi ) { offspring.population[2*i+1].Inversion(rnd); }
          
       }

    Pop.SetPopulation(offspring);
    //Pop.Sort(ProblemSet); //senza MergeSort il codice diventa piuttosto lento
    Pop.MergeSort(Pop, ProblemSet, 0, n_chromosomes - 1);


    }

    //STAMPA popolazione
    void printPopulation() {
        for (int i = 0; i < n_chromosomes; i++) {
            for (int j = 0; j < n_genes; j++) { cout << population[i].GetGenes()[j] << "  "; }
            cout << endl;
        }
    }


    //CROSSOVER
    void Crossover(Random & rnd, int m, int n) { //modifica entrambi i cromosomi

       int i_start = floor(rnd.Rannyu(1.,n_genes)); //indice da cui far partire il crossover
       int support_index=0;

       vector<int> parent1_genes = population[m].GetGenes();
       vector<int> parent2_genes = population[n].GetGenes();

       for(int j=1; j<n_genes; j++) { 
          if( count(parent1_genes.begin()+1, parent1_genes.begin()+i_start, population[n].Get_i_gen(j)) == 0 ) { //controllo che l'elemento j di parent2 non sia nella parte fissata di parent1
             parent1_genes[i_start+support_index]=population[n].Get_i_gen(j); //in tal caso procedo a mettere secondo l'ordine in cui compaiono in parent2
             support_index++;                
          }              
       }

       support_index=0;
       for(int j=1; j<n_genes; j++) {
          if( count(parent2_genes.begin()+1, parent2_genes.begin()+i_start, population[m].Get_i_gen(j)) == 0 ) { 
             parent2_genes[i_start+support_index]=population[m].Get_i_gen(j);
             support_index++;                
          }              
       }

       population[m].SetGenes(parent1_genes);
       population[n].SetGenes(parent2_genes);
       
    }

    

private:
    vector<Chromosome> population;
    int n_chromosomes, n_genes;
    double p; //entra come parametro nella selezione
    double Pc; //Probabilità Crossover
    double Pp; //Probabilità permutazione
    double Ps; //Probabilità shift
    double Ppm; //Probabilità permutazione m città contigue
    double Pi; //Probabilità inversione
};



//FUNZIONI DI STAMPA

//1. Stampa le città nell'ordine in cui sono attraversate
void PrintCities(TSP ProblemSet, Population Pop, int generation, string Problem_type, string Chr_dim) {

   string num_str = to_string(generation);
   ofstream Output("output/"+Problem_type+"/"+Chr_dim+"/Cities/Output.gen"+num_str+".txt");

   int dim = ProblemSet.GetCities().size();

   for(int i=0; i<dim; i++) {
      int index_order = Pop.Get_i_Chrom(0).Get_i_gen(i);
      Output << ProblemSet.GetCities()[index_order].GetX() << setw(20) << ProblemSet.GetCities()[index_order].GetY() << endl;
   }

   Output.close();
}

//2. Stampa il valore di fitness del miglior percorso della generazione presente
void PrintBestFitness(TSP ProblemSet, Population Pop, int generation, string Problem_type, string Chr_dim) {

   ofstream Output("output/"+Problem_type+"/"+Chr_dim+"/Output.BestFitness.txt", ios::app);

   Output << generation << setw(20) << ProblemSet.Fitness1(Pop.Get_i_Chrom(0)) << endl;
   Output.close();
}

//3. Stampa la media del fitness dei cinquanta migliori percorsi della generazione presente
void PrintBestAvFitness(TSP ProblemSet, Population Pop, int generation, string Problem_type, string Chr_dim) {

   ofstream Output("output/"+Problem_type+"/"+Chr_dim+"/Output.BestAvFitness.txt", ios::app);

   int n_best = Pop.GetChromosomes().size();
   double AvFit = 0.;


   for(int i=0; i<n_best; i++) {
      AvFit += ProblemSet.Fitness1(Pop.Get_i_Chrom(i));         
   }

   Output << generation << setw(20) << AvFit/n_best << endl;
   Output.close();
}









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
