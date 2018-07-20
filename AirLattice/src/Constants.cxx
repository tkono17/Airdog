/*
  Constants.cxx
*/
#include "AirLattice/Constants.hxx"

double Constants::g=9.8; // kg m/s^2
double Constants::R=8.31446261815324; // J/(K mol)
double Constants::Cp=30.05; // J/(K mol)

double Constants::MolecularWeight_N2=28.013;
double Constants::MolecularWeight_O2=31.999;
double Constants::MolecularWeight_Air=28.9647;

double cm_to_m() {
  return 1.0E-2;
}
