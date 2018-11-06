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

private:
  ScaleSet();
  static ScaleSet* sInstance;

public:
  static ScaleSet* get();
  ~ScaleSet();

  double setScales(double l0, double v0, double rho, 
		   double theta0, double theta1);

  double toDimLess(double x, ScaleSet::PropertyType tp);
  double fromDimLess(double x, ScaleSet::PropertyType tp);

  double l() const { return mL; }
  double t() const { return mT; }
  double v() const { return mV; }
  double p() const { return mP; }
  double theta0() const { return mTheta0; }
  double theta1() const { return mTheta1; }
  double acceleration() const { return mAcceleration; }

  double kappa() const { return mKappa; }
  double Cp() const { return mCp; }
  double rho() const { return mRho; }
  double beta() const { return mBeta; }
  double g() const { return mG; }
  double Re() const { return mRe; }
  double Pr() const { return mPr; }
  double buoyancyCoefficient() const;

  void print() const;

protected:
  double mL; // Length 
  double mT; // Time
  double mV; // Velocity
  double mRho; // Density
  double mP; // Pressure
  double mTheta0; // Temperature Theta=(theta - theta0)/theta1
  double mTheta1; // Temperature
  double mAcceleration; // Acceleration
  // Parameters
  double mKappa; // Heat conductivity
  double mCp; // Heat capacity / kg
  double mNu; // Diffusion coeffient
  double mBeta; // Coefficient of volume expansion
  double mG; // Gravitational acceleration g
  double mRe; // Reynolds number
  double mPr; // Prandl number
};

#endif // __ScaleSet_hxx__
