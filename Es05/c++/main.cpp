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
#include <vector>
#include <cmath>
#include <functional> //serve per definire funzione di funzione
#include <iomanip> //serve per setprecision
#include "random.h"
#include "stat.h"


using namespace std;

const double PI = 3.14159265358979323846;
const double a0 = 1; // Bohr radius

//FUNZIONI
double Ground(double r, double theta, double phi){
    double psi = pow(M_E,-r)/sqrt(M_PI);
    return psi*psi;
}

double Excited(double r, double theta, double phi){
    double psi = r*pow(M_E,-r/2)*cos(theta)*sqrt(2/M_PI)/8.;
    return psi*psi;
}

//OTHER FUNCTIONS

double min (double a, double b) {
   if (a<b) {return a;}
   else {return b;}
}

double Norm(vector<double> v) {
   int n = v.size(); //dimensione del vettore
   double var = 0.;
   for (int i=0; i<n; i++) {
      var+=v[i]*v[i];
   }
   return sqrt(var);
}

//Trasforma le coordinate cartesiane in coordinate sferiche
vector<double> OrtToSpher(double x, double y, double z) {
    vector<double> SphericalCoord(3);
    SphericalCoord[0] = sqrt(x*x + y*y + z*z );
    SphericalCoord[1] = acos( z/sqrt(x*x + y*y + z*z ) );
    SphericalCoord[2] = (signbit(y) ? -1 : 1) * acos(x/sqrt(x*x + y*y)); //il primo termine restituisce il segno di y[1]
    return SphericalCoord;
}
//Stessa funzione con diversa forma
vector<double> OrtToSpher(vector<double> x) {
    vector<double> SphericalCoord(3);
    SphericalCoord[0] = sqrt(x[0]*x[0] + x[1]*x[1] +x[2]*x[2] );
    SphericalCoord[1] = acos( x[2]/sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2] ) );
    SphericalCoord[2] = (signbit(x[1]) ? -1 : 1) * acos(x[0]/sqrt(x[0]*x[0] + x[1]*x[1])); //il primo termine restituisce il segno di y[1]
    return SphericalCoord;
}


/*FONDAMENTALE: la variabile rnd deve essere passata per riferimento, altrimenti dopo che ha generato numero casuale non vengono aggiornati i parametri
interni e rigenera sempre lo stesso numero casuale!
*/
bool Metropolis(vector<double> &x, double step_size, /*function<double(double,double,double)> p,*/ Random &rnd, bool step, bool groundState) { 
    vector<double> x_Sph = OrtToSpher(x); //Trasformo il punto di partenza in coordinate sferiche
    vector<double> x_proposed(3);
    bool A = false;

    for (int j = 0; j < 3; j++) {
            if (step) { x_proposed[j] = rnd.Gauss(x[j], step_size); }
            else { x_proposed[j] = rnd.Rannyu(x[j]-step_size, x[j]+step_size) ;}
    }

    vector<double> x_proposed_Sph = OrtToSpher(x_proposed);
    double alpha;

    if (!groundState) { alpha = Ground(x_proposed_Sph[0],x_proposed_Sph[1],x_proposed_Sph[2]) / Ground(x_Sph[0],x_Sph[1],x_Sph[2]); } //se 0 (falso) considera ground
    else { alpha = Excited(x_proposed_Sph[0],x_proposed_Sph[1],x_proposed_Sph[2]) / Excited(x_Sph[0],x_Sph[1],x_Sph[2]); }

    if (rnd.Rannyu() < alpha)
    {
        x = x_proposed;
        A = true;
    }
  
    return A;
}

void InitializeRand(Random& rnd); //Implementata sotto







int main (int argc, char *argv[]){


   Random rnd; /*dichiari la variabile rnd appartenente alla classe random, presa dai file random.cpp e random.h */
   InitializeRand(rnd);
 
   int State;
   string sState;
   cout << "Do you want to sample the Ground state or the Excited one?" << endl;
   cout << "Insert 0 for the Ground or 1 for the Excited" << endl;
   cin >> State;

   if(State) { sState = "Excited"; } //Utili per la denominazione dei file di output
   else { sState = "Ground"; }


   //ESERCIZIO 5.1

   int N_samples = 1e6;
   int A_Un = 0; //Conta il tasso di accettazione
   int A_Gauss = 0; 
   vector<double> x_Un = {-1.5,2.5,2.}; //Posizione casuale di partenza
   vector<double> x_Gauss = {-1.5,2.5,2.}; 
   vector<double> r_Un(N_samples); //Conserva le distanze radiali per datablocking
   vector<double> r_Gauss(N_samples); 

   vector<double> stepz(2); //stabilisco i parametri di salto in base allo stato da campionare
   if(!State) { stepz[0]=1.4; stepz[1]=0.8;} //la prima componente è associata al campionamento uniforme, la seconda a quello gaussiano
   else { stepz[0]=2.5; stepz[1]=1.7;}
  

   //EQUILIBRAZIONE 
   int n_eq = 400;
   for (int i=0; i<n_eq; i++){ //Ciclo affinchè i numeri siano generati con una distribuzione che si avvicina a quella asintotica

      ofstream outfileEq("output/outfile"+sState+"Eq.txt", ios::app);

      if ( Metropolis(x_Un, stepz[0], rnd, false, State) ) {A_Un++;}
      if ( Metropolis(x_Gauss, stepz[1], rnd, true, State) )  {A_Gauss++;}

      outfileEq << fixed << setprecision(5) << Norm(x_Un) << setw(20) << Norm(x_Gauss) << endl;
   }

   cout << "Acceptance rate (Uniform): " << double(A_Un)/n_eq << endl;
   cout << "Acceptance rate (Gaussian): " << double(A_Gauss)/n_eq << endl;


   //INIZIO SIMULAZIONE
   A_Un = 0; A_Gauss=0;

   for (int i=0; i<N_samples; i++){

      if ( Metropolis(x_Un, stepz[0], rnd, false, State) ) {A_Un++;}
      if ( Metropolis(x_Gauss, stepz[1], rnd, true, State) ) {A_Gauss++;}

      if( i%100 == 0 ) {//Ogni 100 posizioni generate ne stampo una 

         ofstream outfilePosUn("output/outfile"+sState+"PosUn.txt", ios::app);
         ofstream outfilePosGauss("output/outfile"+sState+"PosGauss.txt", ios::app);
  
         outfilePosUn << x_Un[0] << setw(20) << x_Un[1] << setw(20) << x_Un[2] << endl;
         outfilePosGauss << x_Gauss[0] << setw(20) << x_Gauss[1] << setw(20) << x_Gauss[2] << endl;

         outfilePosUn.close();
         outfilePosGauss.close();
      }

      r_Un[i] = Norm(x_Un);
      r_Gauss[i] = Norm(x_Gauss);


   }

   cout << "Acceptance rate (Uniform): " << double(A_Un)/N_samples << endl;
   cout << "Acceptance rate (Gaussian): " << double(A_Gauss)/N_samples << endl;

   
   //tengo salvati nel main solo i vettori contenenti i dati che poi verranno salvati su file. Gli altri creati e distrutti all'interno della funzione MeanAndErr
   int N_int=200;
   int L=N_samples/N_int;	

   vector<double> sum_prog_Un(N_int);
   vector<double> err_prog_Un(N_int);
   vector<double> sum_prog_Gauss(N_int);
   vector<double> err_prog_Gauss(N_int);

   //Calcolo la media a blocchi
   Stat st;
   st.MeanAndErr(r_Un, sum_prog_Un, err_prog_Un, N_int, L);
   st.MeanAndErr(r_Gauss, sum_prog_Gauss, err_prog_Gauss, N_int, L);

   ofstream outfileRUn("output/outfile"+sState+"RUn.txt");
   ofstream outfileRGauss("output/outfile"+sState+"RGauss.txt");

   for (int i = 0; i < N_int; ++i) {
       outfileRUn << sum_prog_Un[i] << setw(20) << err_prog_Un[i] << endl;
       outfileRGauss << sum_prog_Gauss[i] << setw(20) << err_prog_Gauss[i] << endl;
   }

   outfileRUn.close();
   outfileRGauss.close();

  
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
