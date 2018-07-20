/*
  Atmosphere.cxx
*/
#include "AirLattice/Atmosphere.hxx"
#include "AirLattice/AirProperty.hxx"

Atmosphere::Atmosphere() {
  mFieldData = 0;
}

Atmosphere::~Atmosphere() {
  freeData();
}

int Atmosphere::allocateData(int nx, int ny, int nz) {
  if (mFieldData) {
    freeData();
  }
  mNPoints[0] = nx;
  mNPoints[1] = ny;
  mNPoints[2] = nz;

  int ix, iy;

  mFieldData = new AirProperty**[mNPoints[0]];

  for (ix=0; ix<mNPoints[0]; ++ix) {
    mFieldData[ix] = new AirProperty*[mNPoints[1]];
    for (iy=0; iy<mNPoints[1]; ++iy) {
      mFieldData[ix][iy] = new AirProperty[mNPoints[2]];
    }
  }
  return 0;
}

int Atmosphere::freeData() {
  int ix, iy;

  if (mFieldData) {
    for (ix=0; ix<mNPoints[0]; ++ix) {
      if (mFieldData[ix] == 0) continue;
      for (iy=0; iy<mNPoints[1]; ++iy) {
	if (mFieldData[ix][iy]) {
	  delete [] mFieldData[ix][iy];
	  mFieldData[ix][iy] = 0;
	}
      }
      delete [] mFieldData[ix];
      mFieldData[ix] = 0;
    }
    delete [] mFieldData;
  }
  return 0;
}

AirProperty& Atmosphere::propertyAt(int ix, int iy, int iz) {
  return mFieldData[ix][iy][iz];
}

const AirProperty& Atmosphere::propertyAt(int ix, int iy, int iz) const {
  return mFieldData[ix][iy][iz];
}

void Atmosphere::dump() {
}


