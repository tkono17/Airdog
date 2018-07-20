#ifndef __SystemControl_hxx__
#define __SystemControl_hxx__
/*
  SystemControl.hxx
*/
#include "AirLattice/Environment.hxx"
#include "AirLattice/EnvRecorder.hxx"

class SystemControl {
public:
  SystemControl();
  SystemControl(Environment* env, double dt, double ntimepoints);
  ~SystemControl();

  void setNTimePoints(int x) { mNTimePoints = x; }

  void setDeltaTime(double x) { mDeltaTime = x; }

  void setRecorder(EnvRecorder* rec) { mRecorder = rec; }

  void setRecordInterval(int n) { mRecordInterval = n; }

  void setEnvironment(Environment* env) { mEnvironment = env; }
  Environment* environment() { return mEnvironment; }
  const Environment* environment() const { return mEnvironment; }

  void updateEnvironment(float dt);

  int initialize();
  int run();

protected:
  void updateDensity(); // Continuity equation
  void updateTemperature(); // Heat flow
  void updatePressure(); // State equation
  void updateVelocity(); // Navier-Stokes equation

private:
  double mNTimePoints;
  double mDeltaTime;
  int mTimeStep;
  int mRecordInterval;

  Environment* mEnvironment;
  EnvRecorder* mRecorder;

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
};

#endif // __SystemControl_hxx__
