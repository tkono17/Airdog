/*
  Environment.cxx
*/
#include <iostream>
#include "AirLattice/Environment.hxx"

Environment::Environment() {
  for (int i=0; i<3; ++i) {
    mNPoints[i] = 0;
    mSystemSize[i] = 0.0;
    mElementSize[i] = 0.0;
  }
  mTimeStep = -1;
}

Environment::~Environment() {
  freeData();
}

int Environment::freeData() {
  mAtmosphere.freeData();
  mSurfaceState.freeData();
  return 0;
}

void Environment::setSystemSize(double sx, double sy, double sz) {
  mSystemSize[0] = sx;
  mSystemSize[1] = sy;
  mSystemSize[2] = sz;
}

void Environment::setNPoints(int nx, int ny, int nz) {
  mNPoints[0] = nx;
  mNPoints[1] = ny;
  mNPoints[2] = nz;
}

int Environment::initialize() {
  int i;

  for (i=0; i<3; ++i) {
    mElementSize[i] = mSystemSize[i]/mNPoints[i];
  }

  std::cout << "Initializing the environment:" << std::endl;
  std::cout << "  System size (m): ("
	    << mSystemSize[0] << ", " 
	    << mSystemSize[1] << ", " 
	    << mSystemSize[2] << ")" << std::endl;
  std::cout << "  Lattice structure: ("
	    << mNPoints[0] << ", " 
	    << mNPoints[1] << ", " 
	    << mNPoints[2] << ")" << std::endl;
  std::cout << "  Lattice element size (m): ("
	    << mElementSize[0] << ", " 
	    << mElementSize[1] << ", " 
	    << mElementSize[2] << ")" << std::endl;

  mAtmosphere.allocateData(mNPoints[0], mNPoints[1], mNPoints[2]);
  mSurfaceState.allocateData(mNPoints[0], mNPoints[1]);

  return 0;
}

double Environment::propertyValue(int ix, int iy, int iz, 
				  PropertyType pt) const {
  double value = -1.0;
  if (pt < kSurf_T) {
    const AirProperty& p1 = mAtmosphere.propertyAt(ix, iy, iz);
    if (pt == kRho) {
      value = p1.rho();
    } else if (pt == kP) {
      value = p1.P();
    } else if (pt == kT) {
      value = p1.T();
    } else if (pt == kU) {
      value = p1.u();
    } else if (pt == kV) {
      value = p1.v();
    } else if (pt == kW) {
      value = p1.w();
    } else if (pt == kVaporPressure) {
      value = p1.vaporPressure();
    } else if (pt == kWaterDensity) {
      value = p1.waterDensity();
    }
  } else if (pt < kUnknown) {
    const SurfaceProperty& p2 = mSurfaceState.propertyAt(ix, iy);
    if (pt == kSurf_T) {
      value = p2.T();
    }
  }
  return value;
}

double Environment::rho(int ix, int iy, int iz) const {
  double value=0.0;
  return value;
}

double Environment::P(int ix, int iy, int iz) const {
  double value=0.0;
  ix = (ix + mNPoints[0])%mNPoints[0];
  iy = (iy + mNPoints[1])%mNPoints[1];
  if (iz >= 0 && iz <= (mNPoints[2]-1) ) {
    const AirProperty& p1 = mAtmosphere.propertyAt(ix, iy, iz);
    value = p1.P();
  } else if (iz==-1) {
    const AirProperty& p1 = mAtmosphere.propertyAt(ix, iy, 0);
    value = p1.P();
  } else if (iz == mNPoints[2]) {
    const AirProperty& p1 = mAtmosphere.propertyAt(ix, iy, mNPoints[2]-1);
    value = p1.P();
  }
  return value;
}

double Environment::theta(int ix, int iy, int iz) const {
  double value=0.0;
  ix = (ix + mNPoints[0])%mNPoints[0];
  iy = (iy + mNPoints[1])%mNPoints[1];
  if (iz >= 0 && iz <= (mNPoints[2]-1) ) {
    const AirProperty& p1 = mAtmosphere.propertyAt(ix, iy, iz);
    value = p1.T();
  } else if (iz==-1) {
    const AirProperty& p1 = mAtmosphere.propertyAt(ix, iy, 0);
    const SurfaceProperty& p2 = mSurfaceState.propertyAt(ix, iy);
    value = 2*p2.T()-p1.T();
  } else if (iz == mNPoints[2]) {
    const AirProperty& p1 = mAtmosphere.propertyAt(ix, iy, mNPoints[2]-1);
    value = p1.T();
  }
  return value;
}

double Environment::u(int ix, int iy, int iz) const {
  double value=0.0;
  ix = (ix + mNPoints[0])%mNPoints[0];
  iy = (iy + mNPoints[1])%mNPoints[1];
  if (iz >= 0 && iz <= (mNPoints[2]-1) ) {
    const AirProperty& p1 = mAtmosphere.propertyAt(ix, iy, iz);
    value = p1.u();
  } else if (iz==-1) {
    const AirProperty& p1 = mAtmosphere.propertyAt(ix, iy, 0);
    value = -p1.u();
  } else if (iz == mNPoints[2]) {
    const AirProperty& p1 = mAtmosphere.propertyAt(ix, iy, mNPoints[2]-1);
    value = -p1.u();
  }
  return value;
}

double Environment::v(int ix, int iy, int iz) const {
  double value=0.0;
  ix = (ix + mNPoints[0])%mNPoints[0];
  iy = (iy + mNPoints[1])%mNPoints[1];
  if (iz >= 0 && iz <= (mNPoints[2]-1) ) {
    const AirProperty& p1 = mAtmosphere.propertyAt(ix, iy, iz);
    value = p1.v();
  } else if (iz==-1) {
    const AirProperty& p1 = mAtmosphere.propertyAt(ix, iy, 0);
    value = -p1.v();
  } else if (iz == mNPoints[2]) {
    const AirProperty& p1 = mAtmosphere.propertyAt(ix, iy, mNPoints[2]-1);
    value = -p1.v();
  }
  return value;
}

double Environment::w(int ix, int iy, int iz) const {
  double value=0.0;
  ix = (ix + mNPoints[0])%mNPoints[0];
  iy = (iy + mNPoints[1])%mNPoints[1];
  if (iz >= 0 && iz < (mNPoints[2]-1) ) {
    const AirProperty& p1 = mAtmosphere.propertyAt(ix, iy, iz);
    value = p1.w();
  } else if (iz==-1) {
    value = 0.0;
  } else if (iz >= (mNPoints[2]-1) ) {
    value = 0.0;
  }
  return value;
}

void Environment::dump() {
}


