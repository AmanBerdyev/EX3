#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "ConfigFile.tpp" // Fichier .tpp car inclut un template

using namespace std;



class Exercice3
{

private:
  double t, dt, tFin;
  double m, g, L;
  double d, Omega, kappa;
  double theta, thetadot, thetaOld,thetadotOld;
  int sampling;
  int last;
  ofstream *outputFile;

  void printOut(bool force)
  {
    if((!force && last>=sampling) || (force && last!=1))
    {
      double emec = m *pow(L*thetadot,2) / 2 - m*g*L*(cos(theta));                              // TODO: Evaluer l'energie mecanique
      double pnc  = L*thetadot *(m*d*pow(Omega,2)*sin(Omega * t)-kappa*L*thetadot)*sin(theta);  // TODO: Evaluer la puissance des forces non conservatives

      *outputFile << t << " " << theta << " " << thetadot << " " << emec << " " << pnc << endl;
      last = 1;
    }
    else
    {
      last++;
    }
  }

  double a( double x, double v, double t){
    return (d*pow(Omega,2)*sin(Omega*t)-g)*sin(x)/L - kappa/m*v;
  }

  void step()
  {
    thetaOld = theta;
    thetadotOld = thetadot;

    double A(a(theta,thetadot,t));
    //cout << "a("<< theta << ", " << thetadot << ", " << t << ") = " << A<<endl;

    theta    += thetadot * dt + A * pow(dt,2) /2.0;
    thetadot += A*dt/2.0;
    thetadot  = thetadotOld + (a(thetaOld,thetadot,t) +a(theta,thetadot,t+dt))*dt/2.0;

  }


public:

  Exercice3(int argc, char* argv[])
  {
    string inputPath("configuration.in"); // Fichier d'input par defaut
    if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice3 config_perso.in")
      inputPath = argv[1];

    ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

    for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice3 config_perso.in input_scan=[valeur]")
      configFile.process(argv[i]);

    tFin     = configFile.get<double>("tFin");
    dt       = configFile.get<double>("dt");
    d        = configFile.get<double>("d");
    Omega    = configFile.get<double>("Omega");
    kappa    = configFile.get<double>("kappa");
    m        = configFile.get<double>("m");
    g        = configFile.get<double>("g");
    L        = configFile.get<double>("L");
    theta    = configFile.get<double>("theta0");
    thetadot = configFile.get<double>("thetadot0");
    sampling = configFile.get<int>("sampling");

    // Ouverture du fichier de sortie
    outputFile = new ofstream(configFile.get<string>("output").c_str());
    outputFile->precision(15);
  };

  ~Exercice3()
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
  Exercice3 engine(argc, argv);
  engine.run();
  return 0;
}



