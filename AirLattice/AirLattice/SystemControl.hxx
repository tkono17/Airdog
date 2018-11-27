#ifndef __SystemControl_hxx__
#define __SystemControl_hxx__
/*
  SystemControl.hxx
*/
#include <set>
#include <string>
#include <functional>
#include "AirLattice/Environment.hxx"
#include "AirLattice/EnvRecorder.hxx"

class SystemControl {
public:
  enum boundary_t {
    kXl, kXh, kYl, kYh, kZl, kZh
  };

public:
  SystemControl();
  SystemControl(Environment* env, double dt, double ntimepoints);
  ~SystemControl();

  double elementVolume() const { return mDeltaX*mDeltaY*mDeltaZ; }

  void setNTimePoints(int x) { mNTimePoints = x; }

  void setDeltaTime(double x) { mDeltaTime = x; }

  void setRecorder(EnvRecorder* rec) { mRecorder = rec; }

  void setRecordInterval(int n) { mRecordInterval = n; }

  void setEnvironment(Environment* env) { mEnvironment = env; }
  Environment* environment() { return mEnvironment; }
  const Environment* environment() const { return mEnvironment; }

  void setInitialConditions(const std::set<std::string>& conds) {
    mInitialConditions = conds;
  }
  void setBoundaryConditions(const std::set<std::string>& conds) {
    mBoundaryConditions = conds;
  }
  const std::set<std::string>& initialConditions() const {
    return mInitialConditions;
  }
  const std::set<std::string>& boundaryConditions() const {
    return mBoundaryConditions;
  }

  virtual void updateEnvironment(float dt);

  virtual int initialize();
  virtual int run();

  bool setDoPrint(bool x) { mDoPrint = x; }

protected:
  virtual void applyInitialConditions();
  virtual void updateDensity(); // Continuity equation
  virtual void updateTemperature(); // Heat flow
  virtual void updatePressure(); // Equation of state
  virtual void updateVelocity(); // Navier-Stokes equation
  virtual void applyBoundaryConditions();
  virtual void updateVariables();

  virtual void loop(std::mem_fun_t<void, SystemControl> action);

  void uniformHotspot1();

  double d(std::const_mem_fun_t<double, AirProperty> var, 
	   AirProperty* s0, AirProperty* sm, AirProperty* sp, 
	   double dx, 
	   double mBoundaryValue=0.0, double pBoundaryValue=0.0);

  double dForward(std::const_mem_fun_t<double, AirProperty> var, 
		  AirProperty* s0, AirProperty* sm, AirProperty* sp, 
		  double dx, 
		  double mBoundaryValue=0.0, double pBoundaryValue=0.0);
  
  double dBackward(std::const_mem_fun_t<double, AirProperty> var, 
		   AirProperty* s0, AirProperty* sm, AirProperty* sp, 
		   double dx, 
		   double mBoundaryValue=0.0, double pBoundaryValue=0.0);

  double d2(std::const_mem_fun_t<double, AirProperty> var, 
	    AirProperty* s0, AirProperty* sm, AirProperty* sp, 
	    double dx, 
	    double mBoundaryValue=0.0, double pBoundaryValue=0.0);

  double vd(std::const_mem_fun_t<double, AirProperty> var, 
	    AirProperty* s0, AirProperty* sm, AirProperty* sp, 
	    double v, double dx, 
	    double mBoundaryValue=0.0, double pBoundaryValue=0.0);

  double vd(std::const_mem_fun_t<double, AirProperty> var, double bzl, double bzh);

  double laplacian(std::const_mem_fun_t<double, AirProperty> var, 
		   double bzl, double bzh);

  double d_d(std::const_mem_fun_t<double, AirProperty> var, 
	     AirProperty* s0, AirProperty* sm, AirProperty* sp, 
	     double dx, 
	     double mBoundaryValue=0.0, double pBoundaryValue=0.0);

  double vd_d(std::const_mem_fun_t<double, AirProperty> var, 
	      AirProperty* s0, AirProperty* sm, AirProperty* sp, 
	      double v, double dx, 
	      double mBoundaryValue=0.0, double pBoundaryValue=0.0);

  double vd_d(std::const_mem_fun_t<double, AirProperty> var, double bzl, double bzh);

  bool isAtInnerSite() const;
  bool isAtBoundarySite(boundary_t bt) const;

protected:
  int mNTimePoints;
  double mDeltaTime;
  int mTimeStep;
  int mRecordInterval;

  Environment* mEnvironment;
  EnvRecorder* mRecorder;

  std::set<std::string> mInitialConditions;
  std::set<std::string> mBoundaryConditions;

  // Common derived quantities
  double mAreaYZ;
  double mAreaZX;
  double mAreaXY;
  double mDeltaX;
  double mDeltaY;
  double mDeltaZ;
  // Temporary data shared among functions
  int mXYZIndex[3];
  AirProperty* mAP0;
  AirProperty* mAPx1;
  AirProperty* mAPx2;
  AirProperty* mAPy1;
  AirProperty* mAPy2;
  AirProperty* mAPz1;
  AirProperty* mAPz2;

  bool mDoPrint;
};

#endif // __SystemControl_hxx__
