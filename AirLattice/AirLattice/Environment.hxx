#ifndef __Environment_hxx__
#define __Environment_hxx__
/*
  Environment.hxx
*/
#include "AirLattice/AirProperty.hxx"
#include "AirLattice/Atmosphere.hxx"
#include "AirLattice/SurfaceState.hxx"
#include "AirLattice/PropertyType.hxx"

class Environment {
public:
  Environment();
  ~Environment();

  Atmosphere& atmosphere() { return mAtmosphere; }
  SurfaceState& surfaceState() { return mSurfaceState; }

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

  int freeData();

  void dump();

private:
  int mNPoints[3];
  double mSystemSize[3];
  double mElementSize[3];

  Atmosphere mAtmosphere;
  SurfaceState mSurfaceState;
  int mTimeStep;
};

#endif // __Environment_hxx__
