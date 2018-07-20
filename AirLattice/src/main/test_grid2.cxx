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

  env.setSystemSize(100.0, 1.0, 100.0);
  env.setNPoints(1000, 1, 1000);

  sc.setEnvironment(&env);
  sc.setRecorder(&recorder);
  sc.setDeltaTime(60);
  sc.setNTimePoints(1000);
  sc.setRecordInterval(100);

  env.initialize();
  recorder.openFile("envrec.root");

  sc.run();

  recorder.closeFile();

  return 0;
}
