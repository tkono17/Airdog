/*
  BoundaryCell.cxx
*/
#include "AirLattice/BoundaryCell.hxx"
#include "AirLattice/Environment.hxx"

BoundaryCell::BoundaryCell() {
  mEnv = 0;
  mIx = -1;
  mIx = -1;
  mIx = -1;
  mSurface=0;
  mSurfaceProperty = 0;
  mLocationFlags = 0;
}

BoundaryCell::BoundaryCell(Environment* env, int ix, int iy, int iz) {
  mEnv = env;
  mIx = ix;
  mIx = iy;
  mIx = iz;
  BoundaryCell::Location loc = env->surfaceLocation(ix, iy, iz);
  mSurface = env->surfaceState(loc);
  mSurfaceProperty = 0;
  mLocationFlags = 0;
}

BoundaryCell::~BoundaryCell() {
}

double BoundaryCell::P() const {
  double value=0.0;

  if (mSurface) {
    const AirProperty& ap = innerProperty();
    value = ap.P();
  } else {
    // periodic boundary condition if no surface state is specified
    const AirProperty& ap = periodicProperty();
    value = ap.P();
  }
  return value;
}

double BoundaryCell::theta() const {
  double value = 0.0;
  if (mSurface) {
    const AirProperty& ap = innerProperty();
    value = ap.T();
  } else {
    // periodic boundary condition if no surface state is specified
    const AirProperty& ap = periodicProperty();
    value = ap.T();
  }
  return value;
}

double BoundaryCell::u() const {
  double value = 0.0;
  return value;
}

double BoundaryCell::v() const {
  double value = 0.0;
  return value;
}

double BoundaryCell::w() const {
  double value = 0.0;
  return value;
}

void BoundaryCell::setLocationFlags() {
  const int* np = mEnv->NPoints();

  mLocationFlags = 0;
  if (mIx < 0) mLocationFlags |= (1<<kXLow);
  if (mIx >= np[0]) mLocationFlags |= (1<<kXHigh);
  if (mIy < 0) mLocationFlags |= (1<<kYLow);
  if (mIy >= np[1]) mLocationFlags |= (1<<kYHigh);
  if (mIz < 0) mLocationFlags |= (1<<kZLow);
  if (mIz >= np[2]) mLocationFlags |= (1<<kZHigh);
}

const AirProperty& BoundaryCell::innerProperty() const {
  Atmosphere& atmosphere = mEnv->atmosphere();
  if (isXLow() ) return atmosphere.propertyAt(mIx+1, mIy, mIz);
  if (isXHigh() ) return atmosphere.propertyAt(mIx-1, mIy, mIz);
  if (isYLow() ) return atmosphere.propertyAt(mIx, mIy+1, mIz);
  if (isYHigh() ) return atmosphere.propertyAt(mIx, mIy-1, mIz);
  if (isZLow() ) return atmosphere.propertyAt(mIx, mIy, mIz+1);
  if (isZHigh() ) return atmosphere.propertyAt(mIx, mIy, mIz-1);
  return atmosphere.propertyAt(0, 0, 0);;
}

const AirProperty& BoundaryCell::periodicProperty() const {
  Atmosphere& atmosphere = mEnv->atmosphere();
  const int* np = mEnv->NPoints();

  if (isXLow() ) return atmosphere.propertyAt(mIx+np[0], mIy, mIz);
  if (isXHigh() ) return atmosphere.propertyAt(mIx-np[0], mIy, mIz);
  if (isYLow() ) return atmosphere.propertyAt(mIx, mIy+np[1], mIz);
  if (isYHigh() ) return atmosphere.propertyAt(mIx, mIy-np[1], mIz);
  if (isZLow() ) return atmosphere.propertyAt(mIx, mIy, mIz+np[2]);
  if (isZHigh() ) return atmosphere.propertyAt(mIx, mIy, mIz-np[2]);
  return atmosphere.propertyAt(0, 0, 0);;
}

