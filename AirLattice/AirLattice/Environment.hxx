#ifndef __Environment_hxx__
#define __Environment_hxx__
/*
  Environment.hxx
*/
#include "AirLattice/AirProperty.hxx"
#include "AirLattice/Atmosphere.hxx"
#include "AirLattice/SurfaceState.hxx"
#include "AirLattice/PropertyType.hxx"
#include "AirLattice/BoundaryCondition.hxx"
#include "AirLattice/BoundaryCell.hxx"

class Environment {
public:
  Environment();
  ~Environment();

  Atmosphere& atmosphere() { return mAtmosphere; }
  SurfaceState& surfaceState() { return *(mSurfaceStates[BoundaryCell::kZLow]); }
  SurfaceState* surfaceState(BoundaryCell::Location loc) {
    return mSurfaceStates[loc];
  }

  void setNPoints(int nx, int ny, int nz);
  void setSystemSize(double sx, double sy, double sz);

  double x(int index) const { return mElementSize[0]*index; }
  double y(int index) const { return mElementSize[1]*index; }
  double z(int index) const { return mElementSize[2]*index; }

  const int* NPoints() const { return mNPoints; }
  const double* systemSize() const { return mSystemSize; }
  const double* elementSize() const { return mElementSize; }

  int initialize();

  double propertyValue(int ix, int iy, int iz, PropertyType pt) const;
  double rho(int ix, int iy, int iz) const;
  double P(int ix, int iy, int iz) const;
  double theta(int ix, int iy, int iz) const;
  double u(int ix, int iy, int iz) const;
  double v(int ix, int iy, int iz) const;
  double w(int ix, int iy, int iz) const;

  bool isBoundary(int ix, int iy, int iz) const;

  BoundaryCell::Location surfaceLocation(int ix, int iy, int iz) const;

  const SurfaceProperty& surfacePropertyAt(int ix, int iy, int iz) const;

  const AirProperty& innerPropertyAt(int ix, int iy, int iz) const;

  int freeData();

  void dump();

private:
  int mNPoints[3];
  double mSystemSize[3];
  double mElementSize[3];

  Atmosphere mAtmosphere;
  SurfaceState* mSurfaceStates[6];
  int mTimeStep;
  BoundaryCondition mBoundaryCondition; // Boundary condition flags
};

#endif // __Environment_hxx__
