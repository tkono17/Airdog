/*
  SystemControlMAC.cxx
*/
#include <iostream>
#include "AirLattice/SystemControl.hxx"
#include "AirLattice/SystemControlMAC.hxx"

SystemControlMAC::SystemControlMAC() {
}

SystemControlMAC::SystemControlMAC(Environment* env, double dt, double ntimepoints) {
}

SystemControlMAC::~SystemControlMAC() {
}

void SystemControlMAC::updateEnvironment(float dt) {
  loop(std::mem_fun(&SystemControlMAC::updatePressure) );
  loop(std::mem_fun(&SystemControlMAC::updateVelocity) );
  loop(std::mem_fun(&SystemControlMAC::updateTemperature) );
  //  loop(std::mem_fun(&SystemControl::updateVariables) );

  Atmosphere& atmosphere = mEnvironment->atmosphere();
  std::cout << "Field at (ix, iy, iz)=(49-51, 0, 1)" << std::endl;
  atmosphere.propertyAt(49, 0, 1).print();
  atmosphere.propertyAt(50, 0, 1).print();
  atmosphere.propertyAt(51, 0, 1).print();
  std::cout << "Field at (ix, iy, iz)=(49-51, 0, 0)" << std::endl;
  atmosphere.propertyAt(49, 0, 0).print();
  atmosphere.propertyAt(50, 0, 0).print();
  atmosphere.propertyAt(51, 0, 0).print();
}

int SystemControlMAC::initialize() {
  return SystemControl::initialize();
}

int SystemControlMAC::run() {
  return SystemControl::run();
}

void SystemControlMAC::applyInitialConditions() {
  return uniformHotspot2();
}

void SystemControlMAC::updateTemperature() {
  SystemControl::updateTemperature();
}

void SystemControlMAC::updatePressure() {
  SystemControl::updatePressure();
}

void SystemControlMAC::updateVelocity() {
  SystemControl::updateVelocity();
}

void SystemControlMAC::applyBoundaryConditions() {
  SystemControl::applyBoundaryConditions();
}

void SystemControlMAC::updateVariables() {
  SystemControl::updateVariables();
}

void SystemControlMAC::loop(std::mem_fun_t<void, SystemControlMAC> action) {
  //  SystemControl::loop(action);
}

void SystemControlMAC::uniformHotspot2() {
  uniformHotspot1();
}

