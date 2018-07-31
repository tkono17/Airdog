#ifndef __AirProperty_hxx__
#define __AirProperty_hxx__
/*
  AirProperty.hxx
*/
#include <map>
#include "AirLattice/PropertyType.hxx"

class AirProperty {
public:
  AirProperty();
  ~AirProperty();

  double pressure() const { return mPressure; }
  double P() const { return pressure(); }

  double density() const { return mDensity; }
  double rho() const { return density(); }
  double rhoSI() const { return rho()*1.0E+3; }
  double molDensitySI() const;

  double temperature() const { return mTemperature; }
  double T() const { return temperature(); }

  const double* velocity() const { return mVelocity; }
  double u() const { return mVelocity[0]; }
  double v() const { return mVelocity[1]; }
  double w() const { return mVelocity[2]; }

  double vaporPressure() const { return mVaporPressure; }
  double waterDensity() const { return mWaterDensity; }
  double molecularWeight() const { return mMolecularWeight; }
  double specificHeat() const { return mSpecificHeat; }
  double heatCapacity(double volume) const;
  double heatConductivity() const { return mHeatConductivity; }

  void setPressure(double x) { mPressure = x; }
  void setDensity(double x) { mDensity = x; }
  void setTemperature(double x) { mTemperature = x; }
  void setVelocity(double vx, double vy, double vz);
  void setVelocityU(double vx) { mVelocity[0] = vx; }
  void setVelocityV(double vy) { mVelocity[1] = vy; }
  void setVelocityW(double vz) { mVelocity[2] = vz; }
  void setVaporPressure(double x) { mVaporPressure = x; }
  void setWaterDensity(double x) { mWaterDensity = x; }
  
  void setDensityFromPT();

  // Material characteristic
  void setMolecularWeight(double x) { mMolecularWeight = x; }
  void setSpecificHeat(double x) { mSpecificHeat = x; }
  void setHeatConductivity(double x) { mHeatConductivity = x; }

  void print() const;

  void updateVariable(PropertyType pt, double value);
  void updateVariables();

private:
  // Dynamical properties
  double mPressure; // hPa
  double mDensity; // g/cm3
  double mTemperature; // K
  double mVelocity[3]; // m/s
  double mVaporPressure; // hPa
  double mWaterDensity; // g/cm3

  // Material characteristic
  double mMolecularWeight; // g/mol
  double mSpecificHeat; // Cp [J/(K mol)]
  double mHeatConductivity; // J/(K m2 s)

  // Temporary variable holder
  std::map<int, double> mVariableHolder;

};

#endif // __AirProperty_hxx__
