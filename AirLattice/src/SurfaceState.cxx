/*
  SurfaceState.cxx
*/
#include <iostream>
#include "AirLattice/SurfaceState.hxx"

SurfaceState::SurfaceState() {
  mSurfaceData = 0;
}

SurfaceState::~SurfaceState() {
}

int SurfaceState::allocateData(int nx, int ny) {
  if (mSurfaceData) {
    freeData();
  }
  mNPoints[0] = nx;
  mNPoints[1] = ny;
  int ix;

  std::cout << "Allocating memory for surface data" << std::endl;

  mSurfaceData = new SurfaceProperty*[mNPoints[0]];

  for (ix=0; ix<mNPoints[0]; ++ix) {
    mSurfaceData[ix] = new SurfaceProperty[mNPoints[1]];
  }
  return 0;
}

int SurfaceState::freeData() {
  int ix;

  if (mSurfaceData) {
    for (ix=0; ix<mNPoints[0]; ++ix) {
      if (mSurfaceData[ix] == 0) continue;
      delete [] mSurfaceData[ix];
      mSurfaceData[ix] = 0;
    }
    delete [] mSurfaceData;
  }
  mSurfaceData = 0;
  return 0;
}

SurfaceProperty& SurfaceState::propertyAt(int ix, int iy) {
  return mSurfaceData[ix][iy];
}

const SurfaceProperty& SurfaceState::propertyAt(int ix, int iy) const {
  return mSurfaceData[ix][iy];
}

void SurfaceState::dump() {
  std::cout << "SurfaceState nx, ny=(" << mNPoints[0] << ", " 
	    << mNPoints[1] << ")" << std::endl;
}



