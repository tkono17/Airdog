#ifndef __Constants_hxx__
#define __Constants_hxx__
/*
  Constants.hxx
*/

class Constants {
public:
  // static const double g=9.8; // kg m/s^2
  // static const double R=8.31446261815324; // J/(K mol)
  // static const double Cp=30.05; // J/(K mol)

  // static const double MolecularWeight_N2=28.013;
  // static const double MolecularWeight_O2=31.999;
  // static const double MolecularWeight_Air=28.9647;

  static double g;
  static double R;
  static double T0;
  static double Cp_Air;

  static double MolecularWeight_N2;
  static double MolecularWeight_O2;
  static double MolecularWeight_Air;
};

double cm_to_m();

#endif // __Constants_hxx__
