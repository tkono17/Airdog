#ifndef __SystemControlMAC_hxx__
#define __SystemControlMAC_hxx__
/*
  SystemControlMAC.hxx
*/
#include <set>
#include <string>
#include <vector>
#include <functional>
#include "AirLattice/Environment.hxx"
#include "AirLattice/EnvRecorder.hxx"
#include "AirLattice/SystemControl.hxx"
#include "AirLattice/SparseMatrix.hxx"
#include "AirLattice/LatticeIndexTool.hxx"

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
  void buildMatrixForP(); // build matrix equation to solve P
  void solveP();

  double convU(int ix, int iy, int iz);
  double convV(int ix, int iy, int iz);
  double convW(int ix, int iy, int iz);
  double convTheta(int ix, int iy, int iz);
  double divV(int ix, int iy, int iz);
  double buoyancy(int ix, int iy, int iz);
  double laplacianU(int ix, int iy, int iz);
  double laplacianV(int ix, int iy, int iz);
  double laplacianW(int ix, int iy, int iz);
  double laplacianTheta(int ix, int iy, int iz);

  void loop(std::mem_fun_t<void, SystemControlMAC> action);

  void uniformHotspot2();

  bool isAtInnerSite() const;
  bool isAtBoundarySite(boundary_t bt) const;

private:
  SparseMatrix mM;
  std::vector<double> mB;
  LatticeIndexTool mLatticeIndexTool;

  AirProperty* mAPx1y1;
  AirProperty* mAPx2y1;
  AirProperty* mAPx1y2;
  AirProperty* mAPx2y2;
  //
  AirProperty* mAPz1x1;
  AirProperty* mAPz2x1;
  AirProperty* mAPz1x2;
  AirProperty* mAPz2x2;
  //
  AirProperty* mAPy1z1;
  AirProperty* mAPy2z1;
  AirProperty* mAPy1z2;
  AirProperty* mAPy2z2;

  // Parameters
  double mParameterA;
  double mParamterKappa;
  double mParameterBeta0;
  double mParameterNu;
  double mParameterG;


};

#endif // __SystemControlMAC_hxx__
