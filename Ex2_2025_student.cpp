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
      static int counter = 0; // Contatore statico per mantenere il valore tra le chiamate

      if((!force && last >= sampling) || (force && last != 1))
      {
          double emec = (1.0 / 24.0) * (thetadot * thetadot) * (m * L * L) - mu * B0 * cos(theta);
          double pnc  = -thetadot * mu * sin(theta) * (B1 * sin(Omega * t )) - kappa * thetadot * thetadot;


              *outputFile << t << " " << theta << " " << thetadot << " " << emec << " " << pnc << endl;

          counter++; // Incrementa il contatore

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
        acc[0] = -(mu/I) * sin(x ) * (B0 + B1 * sin(Omega * t_)); //NON CAPISCO LA FORMULA THETA PUNTO DEL OVERLEAF DOVE L'IMPLEMENTO?
        acc[1] = -(1.0/I) * kappa * v; //Deve dipendere da qualcos'altro no, perche cosi per kappa = 0 rimane sempre nullo:
        return acc;
    }


    void step()
  {
    // TODO: implement the extended Verlet scheme Section 2.7.4

    // cout << "Theta " << theta << " Thetadot " << thetadot << " t " << t << "intervalle de temps donnée par" << dt<< endl;

    valarray<double> a_j = acceleration(theta, thetadot, t); //Otteniamo l'accelerazione a(x_j, v_j, t_j)

    // cout << "a(x_j, v_j, t_j) " << a_j[0] << " " << a_j[1] << endl;



    //Position
    theta                = theta + thetadot*dt + (0.5) * (a_j[0] + a_j[1]) * dt*dt; //Calcolo do x_j+1

    //Vitesse de Demi-pas:
    thetadot_demi_pas = thetadot + (0.5) * (a_j[0] + a_j[1]) * dt;

    //Un pas d'accellaration
    valarray<double> a_new = acceleration(theta, 0, t + dt); //qui otteniamo a(x_j+1, v_j+1, t_j+1)
    valarray<double> a_new_new = acceleration(0, thetadot_demi_pas, t + (0.5)*dt); //qui otteniamo a(x_j+1, v_j+1, t_j+1)

    I = (1.0/12.0) * m * L * L;


    thetadot = thetadot + (0.5) * (a_j[0] + a_new[0]) * dt + (a_new_new[1]) * dt;

}




public:
  Exercice2(int argc, char* argv[])
  {
    const double pi=3.1415926535897932384626433832795028841971e0;
    string inputPath("configuration.in.example"); // Fichier d'input par defaut
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

    double I = m*L*L/12;
    om_0 = sqrt(mu*B0/I);
    Omega = 2*om_0;
    cout << "Value om_0 :"<< om_0 << endl;
    cout << "Value Omega :"<< Omega << endl;

    // TODO: implement the expression for tFin, dt
    T = 2 * pi /Omega;

    sampling  = T;

    if(N_excit>0){
      // TODO Définir le pas temporel et le temp final si N_excit est défini.
      cout << "Tempo totale" << T << endl;
      dt= T/nsteps_per;
      tFin= N_excit*T;
    }else{
      cout << "NON FUNZIONA" << endl;
      tFin = T;
      dt = tFin/nsteps_per;
    }


    cout << "final time is "<<"  "<< tFin << endl;

  }

  void sol_analytique()
  {
    double theta0 = 1e-3;
    double theta_fin = theta0 * cos(om_0*tFin);
    double thetaDot_fin = - theta0 * om_0 * sin(om_0*tFin);
    cout << "Solution analytique position : " << theta_fin << endl;
    cout << "Solution analytique vitesse : " << thetaDot_fin << endl;
    cout << "omega_0 :" << om_0 << endl;
  };


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
