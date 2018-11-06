/*
  test_mac1.cxx
*/
#include <iostream>
#include "AirLattice/Environment.hxx"
#include "AirLattice/EnvRecorder.hxx"
#include "AirLattice/SystemControlMAC.hxx"
#include "AirLattice/ScaleSet.hxx"

int main(int argc, char* argv[]) {
  Environment env;
  EnvRecorder recorder(&env);
  SystemControlMAC sc;

  //  env.setSystemSize(100.0E+3, 1.0, 50.0E+3);
  env.setSystemSize(10.0, 1.0, 10.0);
  env.setNPoints(100, 1, 100);

  sc.setEnvironment(&env);
  sc.setRecorder(&recorder);
  //  sc.setDeltaTime(60);
  sc.setDeltaTime(1.0E-4);
  sc.setNTimePoints(10000);
  sc.setRecordInterval(100);

  sc.setNTimePoints(1000);
  sc.setRecordInterval(1);
  sc.setDoPrint(true);

  double l0 = 10.0;
  double v0 = 1.0E+2;
  double rho0 = 1.1839; // kg/m3
  double nu = 0.2E-4; // m2/s
  double theta0 = 273.15;
  double theta1 = 1.0;

  ScaleSet* ss = ScaleSet::get();
  ss->setScales(l0, v0, rho0, theta0, theta1);

  env.initialize();
  recorder.openFile("envrec.root");

  sc.initialize();

  sc.run();

  recorder.closeFile();

  return 0;
}
