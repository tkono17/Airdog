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

  void updateEnvironment(float dt);

  int initialize();
  int run();

protected:
  void applyInitialConditions();
  void updateDensity(); // Continuity equation
  void updateTemperature(); // Heat flow
  void updatePressure(); // Equation of state
  void updateVelocity(); // Navier-Stokes equation
  void applyBoundaryConditions();
  void updateVariables();

  void loop(std::mem_fun_t<void, SystemControl> action);

  void uniformHotspot1();

private:
  double mNTimePoints;
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
