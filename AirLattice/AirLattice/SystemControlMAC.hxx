#ifndef __SystemControlMAC_hxx__
#define __SystemControlMAC_hxx__
/*
  SystemControlMAC.hxx
*/
#include <set>
#include <string>
#include <functional>
#include "AirLattice/Environment.hxx"
#include "AirLattice/EnvRecorder.hxx"
#include "AirLattice/SystemControl.hxx"

class SystemControlMAC : public SystemControl {
public:
  enum boundary_t {
    kXl, kXh, kYl, kYh, kZl, kZh
  };

public:
  SystemControlMAC();
  SystemControlMAC(Environment* env, double dt, double ntimepoints);
  ~SystemControlMAC();

  void updateEnvironment(float dt);

  int initialize();
  int run();

  bool setDoPrint(bool x) { mDoPrint = x; }

protected:
  void applyInitialConditions();
  void updateTemperature(); // Heat flow
  void updatePressure(); // Equation of state
  void updateVelocity(); // Navier-Stokes equation
  void applyBoundaryConditions();
  void updateVariables();

  void loop(std::mem_fun_t<void, SystemControlMAC> action);

  void uniformHotspot2();

  bool isAtInnerSite() const;
  bool isAtBoundarySite(boundary_t bt) const;

private:
  // Parameters
  double mParameterA;
  double mParameterBeta0;
  double mParameterNu;
  double mParameterG;


};

#endif // __SystemControlMAC_hxx__
