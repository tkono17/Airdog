#ifndef __SurfaceState_hxx__
#define __SurfaceState_hxx__
/*
  SurfaceState.hxx
*/
#include "AirLattice/SurfaceProperty.hxx"

class SurfaceState {
public:
  SurfaceState();
  ~SurfaceState();

  int allocateData(int nx, int ny);
  int freeData();

  SurfaceProperty& propertyAt(int ix, int iy);
  const SurfaceProperty& propertyAt(int ix, int iy) const;

  void dump();

private:
  int mNPoints[2];

  SurfaceProperty** mSurfaceData;
};

#endif // __SurfaceState_hxx__
