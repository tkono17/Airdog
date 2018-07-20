/*
  Environment.cxx
*/
#include "AirLattice/Environment.hxx"

Environment::Environment() {
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

void Environment::dump() {
}


