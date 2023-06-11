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
#include <ostream>
#include <cmath>
#include <iomanip>
#include <string>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main() {


  //DATI EQUILIBRAZIONE  
  for (int i=0; i<3; i++) { //Ciclo su diverse temperature su cui si equilibra
    
    Input(); //Ogni volta riparto da una configurazione non equilibrata
    if (i!=0) { 
      temp=temp+i*0.5; 
      beta=1/temp;
    }
    Equilibrazione(temp);

  }


  //SIMULAZIONE
  Input(); //Resetta le variabili di input

  if(restart){ //Termalizzazione. Scegliamo 1000 in base ai dati dell'equilibrazione.
    for (int i = 0; i<1000; i++) Move(metro);
  }

  for(int iTemp=0; iTemp<=10; iTemp++) {

    if (iTemp!=0) {
      temp = temp+DeltaTemp;
      beta = 1./temp;
    }

    for(int iblk=1; iblk <= nblk; ++iblk) {
      Reset(iblk);   //Reset block averages
      for(int istep=1; istep <= nstep; ++istep) {
        Move(metro);
        Measure();
        Accumulate(); //Update block averages
      }
      Averages(iblk);   //Print results for current block
    }

    Averages_temp(temp);

  }

  ConfFinal(); //Write final configuration


  return 0;

}


void Input(void)
{
  ifstream ReadInput, Seed;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();
  
//Read input informations
  ReadInput.open("input.dat");
  ReadInput >> restart; //Se restart è 0 (-->Bool=false) allora si riparte dagli stessi numeri di input. Si rifà simulazione identica
  ReadInput >> startingpoint; 

  if(!restart) Seed.open("seed.out");
  else Seed.open("seed.in");
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed,p1,p2);
  Seed.close();

  ReadInput >> temp;
  beta = 1.0/temp;
  Temperature = to_string(round(temp*10)/10.);
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs
  if (metro) { Method = "Met"; }
  else { Method = "Gibbs"; }

  ReadInput >> nblk;

  ReadInput >> nstep;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

//initial configuration

  if(!startingpoint) { //questo primo blocco if/else potrebbe essere posto sotto if(restart), tuttavia quando si fa la prima simulazione, in caso di restart=0 non c'è nessun file config.final da aprire.
    for (int i=0; i<nspin; ++i) {
      s[i] = 1; //parte da stato ordinato a bassa temperatura
    } 
  }
  else {
    for (int i=0; i<nspin; ++i) {
      if(rnd.Rannyu() >= 0.5) {s[i] = 1;}
      else {s[i] = -1;} //parte da stato disordinato a alta temperatura
    }
  }

  if(!restart) {
    ifstream ReadConfig;
    ReadConfig.open("config.final");
    for (int i=0; i<nspin; ++i) {
      ReadConfig >> s[i]; //prende gli spin dalla configurazione già equilibrata
    } 
  } 
  
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl << endl;
}


void Move(int metro) {

  int o;

  for(int i=0; i<nspin; ++i) {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    if(metro==1) {//Metropolis
       double DeltaE = Boltzmann(-s[o],o)-Boltzmann(s[o],o);
       //double DeltaE = 2*J*s[o]*(s[Pbc(o-1)]+s[Pbc(o+1)]) + 2*h*s[o];
       double p = min( 1. , exp(-beta*(DeltaE)) );
       if (rnd.Rannyu() < p) { //genera un numero casuale anche quando DeltaE è negativo (in tal caso si entra nell'if in modo scontato)
          s[o] = -s[o];
          accepted++;
       }
       attempted++;  
    }

    else {//Gibbs sampling
       s[o]=1; //si potrebbero tenere in considerazione anche spin negativi ma non è necessario (bisognerebbe cambiare segno nella distrib p).         
       double DeltaE = 2.*J*s[o]*(s[Pbc(o-1)]+s[Pbc(o+1)]) + 2*h*s[o];
       double p = 1./(1+exp(-beta*DeltaE));
            
       if (rnd.Rannyu() > p) {s[o] = -1;}
       accepted++;
       attempted++;  
    }
  }
}

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm; //non si potrebbe implementare mettendo s[ip] al posto che sm nell'input?
  return ene;
}

void Measure()
{
  double u = 0.0, m = 0.0;

  //cycle over spins
  for (int i=0; i<nspin; ++i) {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]); //calcolo energia totale
     m += s[i]; //calcolo magnetizzazione totale
  }

  walker[iu] = u;
  walker[im] = m; 
  walker[ic] = u*u;
  walker[ix] = m*m;

}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1) {
       for(int i=0; i<n_props; ++i) {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i) {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i) {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0; // per come è costruito il codice alla fine è uguale a nstep (# di configurazioni per blocco)
}


void Averages(int iblk) //Print results for current block
{
    
   ofstream Ene, Heat, Mag, Chi;
   const int wd=20;
    
    cout << "Temperature " << temp << endl;
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
  
    if (h==0.) {  
       //Ene.open("data/output.ene", ios::app);
       //Ene << std::scientific << std::setprecision(9);
       stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy per spin
       glob_av[iu]  += stima_u;
       glob_av2[iu] += stima_u*stima_u;
       err_u=Error(glob_av[iu],glob_av2[iu],iblk);
       //Ene << temp << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
       //Ene.close();

       //Heat.open("data/output.heat", ios::app);
       //Heat << std::scientific << std::setprecision(9);
       stima_c = pow(beta,2) * ( blk_av[ic]/blk_norm - pow(blk_av[iu]/blk_norm,2) )/(double)nspin; //Specific heat
       glob_av[ic]  += stima_c;
       glob_av2[ic] += stima_c*stima_c;
       err_c=Error(glob_av[ic],glob_av2[ic],iblk);
       //Heat << temp << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
       //Heat.close();

       //Chi.open("data/output.chi", ios::app);
       //Chi << std::scientific << std::setprecision(9);
       stima_x = beta * ( blk_av[ix]/blk_norm - pow(blk_av[im]/blk_norm,2) )/(double)nspin; //Susceptibility
       glob_av[ix]  += stima_x;
       glob_av2[ix] += stima_x*stima_x;
       err_x=Error(glob_av[ix],glob_av2[ix],iblk);
       //Chi << temp << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
       //Chi.close();
    }    

    else {
       //Mag.open("data/output.mag", ios::app);
       //Mag << std::scientific << std::setprecision(9);
       stima_m = blk_av[im]/blk_norm/(double)nspin; //Magnetization per spin
       glob_av[im]  += stima_m;
       glob_av2[im] += stima_m*stima_m;
       err_m=Error(glob_av[im],glob_av2[im],iblk);
       //Mag << temp << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
       //Mag.close();
    }

    cout << "----------------------------" << endl << endl;
}


void Averages_temp(double temp) {
    
   ofstream Ene_t, Heat_t, Mag_t, Chi_t;
   const int wd=20;

   if (metro==1) { 
     if(h==0){
       Ene_t.open("data/energy_temp_metropolis.out", ios::app);
       Ene_t << std::scientific << std::setprecision(9);
       Ene_t << setw(wd) << temp <<  setw(wd) << glob_av[iu]/(double)nblk << setw(wd) << err_u << endl;
       Ene_t.close();

       Heat_t.open("data/heat_temp_metropolis.out", ios::app);
       Heat_t << std::scientific << std::setprecision(9);
       Heat_t << setw(wd) << temp <<  setw(wd) << glob_av[ic]/(double)nblk << setw(wd) << err_c << endl;
       Heat_t.close();
    
       Chi_t.open("data/chi_temp_metropolis.out", ios::app);
       Chi_t << std::scientific << std::setprecision(9);
       Chi_t << setw(wd) << temp << setw(wd) << glob_av[ix]/(double)nblk << setw(wd) << err_x << endl;
       Chi_t.close();
     }
     else if (h!=0) {
       Mag_t.open("data/magnetization_temp_metropolis.out", ios::app);
       Mag_t << std::scientific << std::setprecision(9);
       Mag_t << setw(wd) << temp <<  setw(wd) << glob_av[im]/(double)nblk << setw(wd) << err_m << endl;
       Mag_t.close();
     }
   } 
   else {
     if(h==0){
       Ene_t.open("data/energy_temp_gibbs.out", ios::app);
       Ene_t << std::scientific << std::setprecision(9);
       Ene_t << setw(wd) << temp <<  setw(wd) << glob_av[iu]/(double)nblk << setw(wd) << err_u << endl;
       Ene_t.close();

       Heat_t.open("data/heat_temp_gibbs.out", ios::app);
       Heat_t << std::scientific << std::setprecision(9);
       Heat_t << setw(wd) << temp <<  setw(wd) << glob_av[ic]/(double)nblk << setw(wd) << err_c << endl;
       Heat_t.close();
    
       Chi_t.open("data/chi_temp_gibbs.out", ios::app);
       Chi_t << std::scientific << std::setprecision(9);
       Chi_t << setw(wd) << temp << setw(wd) << glob_av[ix]/(double)nblk << setw(wd) << err_x << endl;
       Chi_t.close();
     } 
     else if(h!=0) {
       Mag_t.open("data/magnetization_temp_gibbs.out", ios::app);
       Mag_t << std::scientific << std::setprecision(9);
       Mag_t << setw(wd) << temp <<  setw(wd) << glob_av[im]/(double)nblk << setw(wd) << err_m << endl;
       Mag_t.close();
     }
   }

}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i) {
    WriteConf << s[i] << endl;
  }

  WriteConf.close();
  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

void Equilibrazione(double temp) {

    Temperature = to_string(round(temp*10)/10.);
    restart = 1;

    ofstream Equilibrazione;
    string Filename_Eq = "data/Equilibrazione" + Method + Temperature;
    Equilibrazione.open(Filename_Eq);
    const int wd=20;

    for(int i=1; i <= 500; ++i) { //Equilibrazione

      Move(metro);
      Measure();

      Equilibrazione << std::scientific << std::setprecision(9);
      stima_u = walker[iu]/(double)nspin; //Energy per spin
      stima_m = walker[im]/(double)nspin; //Magnetization per spin
      Equilibrazione <<  setw(wd) << stima_u << setw(wd) << stima_m << endl;

    }

    Equilibrazione.close();

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
