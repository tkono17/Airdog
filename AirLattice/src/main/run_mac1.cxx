/*
  test_mac1.cxx
*/
#include <iostream>
#include <cstdio>
#include "AirLattice/Environment.hxx"
#include "AirLattice/EnvRecorder.hxx"
#include "AirLattice/SystemControlMAC.hxx"
#include "AirLattice/ScaleSet.hxx"

int main(int argc, char* argv[]) {
  double dT=1.0E-8;
  double NT=1000;

  if (argc > 1) {
    std::sscanf(argv[1], "%e", &dT);
  }
  if (argc > 2) {
    std::sscanf(argv[2], "%e", &NT);
  }

  Environment env;
  EnvRecorder recorder(&env);
  SystemControlMAC sc;

  double l0 = 1.0;
  double v0 = 1.0E+2;
  double rho0 = 1.1839; // kg/m3
  double theta0 = 273.15;
  double theta1 = 1.0;

  ScaleSet* ss = ScaleSet::get();
  //v0 = 1.0E+3;
  v0 = 1.0;
  ss->setScales(l0, v0, rho0, theta0, theta1);
  v0 = ss->nu()/l0;
  ss->setScales(l0, v0, rho0, theta0, theta1);

  //  env.setSystemSize(100.0E+3, 1.0, 50.0E+3);
  env.setSystemSize(1, 0.1, 1);
  env.setNPoints(10, 1, 10);

  sc.setEnvironment(&env);
  sc.setRecorder(&recorder);

  //  sc.setDeltaTime(60);
  sc.setDeltaTime(dT);
  sc.setNTimePoints(NT);
  sc.setRecordInterval(100);

  sc.setNTimePoints(100);
  sc.setRecordInterval(1);
  sc.setDoPrint(true);

  env.initialize();
  recorder.openFile("envrec.root");

  sc.initialize();

  sc.run();

  recorder.closeFile();

  return 0;
}
