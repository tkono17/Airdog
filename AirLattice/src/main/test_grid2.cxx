/*
  test_grid2.cxx
*/
#include <iostream>
#include "AirLattice/Environment.hxx"
#include "AirLattice/EnvRecorder.hxx"
#include "AirLattice/SystemControl.hxx"

int main(int argc, char* argv[]) {
  Environment env;
  EnvRecorder recorder(&env);
  SystemControl sc;

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

  env.initialize();
  recorder.openFile("envrec.root");

  sc.initialize();
  std::cout << "Run" << std::endl;
  sc.run();

  std::cout << "Close recorder file" << std::endl;

  recorder.closeFile();

  return 0;
}
