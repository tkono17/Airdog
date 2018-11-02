#ifndef __ScaleSet_hxx__
#define __ScaleSet_hxx__
/*
  ScaleSet.hxx
*/

class ScaleSet {
public:
  enum PropertyType {
    kL, kT, kV, kP, kParameterA, kTheta, kAcceleration
  };

public:
  ScaleSet();
  ~ScaleSet();

  double setScales(double l0, double v0, double rho, 
		   double theta0, double theta1);

  double toDimLess(double x, ScaleSet::PropertyType tp);
  double fromDimLess(double x, ScaleSet::PropertyType tp);

  double l() const { return mL; }
  double t() const { return mT; }
  double v() const { return mV; }
  double p() const { return mP; }
  double parameterA() const { return mParameterA; }
  double theta0() const { return mTheta0; }
  double theta1() const { return mTheta1; }
  double acceleration() const { return mAcceleration; }

protected:
  double mL; // Length 
  double mT; // Time
  double mV; // Velocity
  double mRho; // Density
  double mP; // Pressure
  double mParameterA; // Parameter a=kappa/(cp*rho)
  double mTheta0; // Temperature Theta=(theta - theta0)/theta1
  double mTheta1; // Temperature
  double mAcceleration; // Acceleration
};

#endif // __ScaleSet_hxx__
