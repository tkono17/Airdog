#ifndef __AirProperty_hxx__
#define __AirProperty_hxx__
/*
  AirProperty.hxx
*/

class AirProperty {
public:
  AirProperty();
  ~AirProperty();

  double pressure() const { return mPressure; }
  double P() const { return pressure(); }

  double density() const { return mDensity; }
  double rho() const { return density(); }

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
  double heatCapacity() const { return mSpecificHeat; }
  double heatConductivity() const { return mHeatConductivity; }

  void setPressure(double x) { mPressure = x; }
  void setDensity(double x) { mDensity = x; }
  void setTemperature(double x) { mTemperature = x; }
  void setVelocity(double vx, double vy, double vz);
  void setVaporPressure(double x) { mVaporPressure = x; }
  void setWaterDensity(double x) { mWaterDensity = x; }
  
  // Material characteristic
  void setMolecularWeight(double x) { mMolecularWeight = x; }
  void setSpecificHeat(double x) { mSpecificHeat = x; }
  void setHeatConductivity(double x) { mHeatConductivity = x; }

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
};

#endif // __AirProperty_hxx__
