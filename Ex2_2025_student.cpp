#include <iostream>
#include <fstream>
#include <iomanip>
#include "ConfigFile.h" // Il contient les methodes pour lire inputs et ecrire outputs 
#include <valarray>
#include <cmath> // Se usi fun
using namespace std;

class Exercice2
{

private:
  double t, dt, tFin;
  double m, g, L;
  double Omega, kappa;
  double theta, thetadot;
  double B0, B1;
  double mu;
  double om_0, om_1, nu, Ig;

  int N_excit, nsteps_per, Nperiod;
  int sampling;
  int last;

  double T;
  double I;


  double thetadot_demi_pas;
  ofstream *outputFile;

  void printOut(bool force)
  {
    if((!force && last>=sampling) || (force && last!=1))
    {
      double emec = (1.0/24.0) * (thetadot * thetadot) *(m * L*L) - mu* B0*cos(theta);
      double pnc  = -thetadot *mu * sin(theta)*( B1*sin(Omega*t)) - kappa*thetadot*thetadot; 

      *outputFile <<  t << " " << theta << " " << thetadot << " " << emec << " " << pnc << endl;
      last = 1;
    }
    else
    {
      last++;
    }
  }

    // TODO define angular acceleration functions and separate contributions
    // for a[0] = function of (x,t), a[1] = function of (v)
    valarray<double> acceleration(double x, double v, double t_)
    {
        valarray<double> acc(2);
        I = (1.0/12.0) * m * L * L;
        acc[0] = -(mu/I) * sin(x) * (B0 + B1 * sin(Omega * t_));
        acc[1] = -(1.0/I) * kappa * v;
        return acc;
    }
  void step() {
      valarray<double> a = acceleration(theta, thetadot, t);
      // TODO Définir les variables dont vous avez besoin et fair le schema numerique

      theta = theta + thetadot * dt + 0.5 * pow(dt, 2) * (a[0] + a[1]);

      while (theta >M_PI or theta<-M_PI) {
          if (theta > M_PI) {
              theta -= 2 * M_PI;
          } else if (theta < -M_PI) {
              theta += 2 * M_PI;
          }
      }



  double thetadot_demi     = thetadot + 0.5*dt*(a[0] + a[1]);

  valarray<double> new_a = acceleration(theta,thetadot_demi, t+dt);
  valarray<double> new_a_demi = acceleration(theta,thetadot_demi, t+0.5*dt);
  thetadot = thetadot + 0.5*(a[0] + new_a[0])*dt + new_a_demi[1]*dt;

}
    


/*     void step()
  {
    // TODO: implement the extended Verlet scheme Section 2.7.4
    
    valarray<double> a = acceleration(theta, thetadot, t);
    cout << "Theta " << theta << " Thetadot " << thetadot << " t " << t << endl; 
    //Position
    theta                = theta + thetadot*dt + (0.5) * (a[0] + a[1]) * dt*dt;
    //Demi-pas:
    thetadot_demi_pas = thetadot + (0.5) * (a[0] + a[1]) * dt;
    cout << "Thetadot_demi_pas " << thetadot_demi_pas << endl;
    //Pas complet
    valarray<double> a_1 = acceleration(theta, thetadot_demi_pas, t + dt);
    cout << "a_1 " << a_1[0] << " " << a_1[1] << endl;
    valarray<double> a_2 = acceleration(theta, thetadot_demi_pas, t + 0.5*dt);
    cout << "a_2 " << a_2[0] << " " << a_2[1] << endl;
    thetadot = thetadot + (0.5) * (a[0] + a_1[0]) * dt + (0.5) * (a[1] + a_1[1]) * dt;
    cout << "Thetadot " << thetadot << endl;

} */




public:
  Exercice2(int argc, char* argv[])
  {
    const double pi=3.1415926535897932384626433832795028841971e0;
    string inputPath("configuration.in"); // Fichier d'input par defaut
    if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice2 config_perso.in")
      inputPath = argv[1];

    ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

    for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice2 config_perso.in input_scan=[valeur]")
      configFile.process(argv[i]);

    Omega    = configFile.get<double>("Omega");     // frequency of oscillating magnetic field
    kappa    = configFile.get<double>("kappa");     // coefficient for friction
    m        = configFile.get<double>("m");         // mass
    L        = configFile.get<double>("L");         // length
    B1       = configFile.get<double>("B1");     //  B1 part of magnetic fields (oscillating part amplitude)
    B0       = configFile.get<double>("B0");     //  B0 part of magnetic fields (static part amplitude)
    mu       = configFile.get<double>("mu");     //  magnetic moment 
    theta    = configFile.get<double>("theta0");    // initial condition in theta
    thetadot = configFile.get<double>("thetadot0"); // initial condition in thetadot
    sampling = configFile.get<int>("sampling");     // number of time steps between two writings on file
    N_excit  = configFile.get<int>("N_excit");      // number of periods of excitation
    Nperiod  = configFile.get<int>("Nperiod");      // number of periods of oscillation of the eigenmode
    nsteps_per= configFile.get<int>("nsteps");      // number of time step per period

    // Ouverture du fichier de sortie
    outputFile = new ofstream(configFile.get<string>("output").c_str());
    outputFile->precision(15);

    // define auxiliary variables if you need/want
    
    // TODO: implement the expression for tFin, dt
    T = 2 * pi /Omega;

    if(N_excit > 0){
      tFin = N_excit * T;
      dt = T / nsteps_per;
    } else if (Nperiod > 0) {
      tFin = Nperiod * 2 * pi / sqrt(mu * B0 / (m * L * L)); // Frequenza naturale del sistema
      dt = tFin / (Nperiod * nsteps_per);
    } else {
      tFin = T;
      dt = T / nsteps_per;
    }
    
    
    cout << "final time is "<<"  "<< tFin << endl; 

  }
  

  ~Exercice2()
  {
    outputFile->close();
    delete outputFile;
  };

    void run()
  {
    t = 0.;
    last = 0;
    printOut(true);

    while( t < tFin-0.5*dt )
    {
      step();
      t += dt;
      printOut(false);
    }
    printOut(true);
  };

};

int main(int argc, char* argv[])
{
  Exercice2 engine(argc, argv);
  engine.run();

  return 0;
}
