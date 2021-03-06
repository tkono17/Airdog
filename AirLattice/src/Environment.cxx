/*
  Environment.cxx
*/
#include <iostream>
#include "AirLattice/Environment.hxx"
#include "AirLattice/BoundaryCell.hxx"

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
  for (int i=0; i<6; ++i) {
    if (mSurfaceStates[i]) {
      mSurfaceStates[i]->freeData();
    }
  }
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
	    << mSystemSize[1]  << ", " 
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

  for (i=0; i<6; ++i) {
    mSurfaceStates[i] = new SurfaceState();
  }
  if (mSurfaceStates[BoundaryCell::kXLow]) {
    mSurfaceStates[BoundaryCell::kXLow]->allocateData(mNPoints[1], mNPoints[2]);
  }
  if (mSurfaceStates[BoundaryCell::kXHigh]) {
    mSurfaceStates[BoundaryCell::kXHigh]->allocateData(mNPoints[1], mNPoints[2]);
  }
  if (mSurfaceStates[BoundaryCell::kYLow]) {
    mSurfaceStates[BoundaryCell::kYLow]->allocateData(mNPoints[2], mNPoints[0]);
  }
  if (mSurfaceStates[BoundaryCell::kYHigh]) {
    mSurfaceStates[BoundaryCell::kYHigh]->allocateData(mNPoints[2], mNPoints[0]);
  }
  if (mSurfaceStates[BoundaryCell::kZLow]) {
    mSurfaceStates[BoundaryCell::kZLow]->allocateData(mNPoints[0], mNPoints[1]);
  }
  if (mSurfaceStates[BoundaryCell::kZHigh]) {
    mSurfaceStates[BoundaryCell::kZHigh]->allocateData(mNPoints[0], mNPoints[1]);
  }

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
    const SurfaceProperty& p2 = mSurfaceStates[BoundaryCell::kZLow]->propertyAt(ix, iy);
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
  if (isBoundary(ix, iy, iz) ) {
    const SurfaceProperty& sp = surfacePropertyAt(ix, iy, iz);
    const AirProperty& ap1 = innerPropertyAt(ix, iy, iz);
    value = ap1.P();
  } else {
    const AirProperty& ap0 = mAtmosphere.propertyAt(ix, iy, iz);
    value = ap0.P();
  }
  // ix = (ix + mNPoints[0])%mNPoints[0];
  // iy = (iy + mNPoints[1])%mNPoints[1];
  // if (iz >= 0 && iz < mNPoints[2]) {
  //   const AirProperty& p1 = mAtmosphere.propertyAt(ix, iy, iz);
  //   value = p1.P();
  // } else if (iz==-1) {
  //   const AirProperty& p1 = mAtmosphere.propertyAt(ix, iy, 0);
  //   value = p1.P();
  // } else if (iz == mNPoints[2]) {
  //   const AirProperty& p1 = mAtmosphere.propertyAt(ix, iy, mNPoints[2]-1);
  //   value = p1.P();
  // }
  return value;
}

double Environment::theta(int ix, int iy, int iz) const {
  double value=0.0;
  // if (ix < 0) {
  //   std::cout << "Get theta at XLow" << std::endl;
  // } else if (ix >= mNPoints[0]) {
  //   std::cout << "Get theta at XHigh" << std::endl;
  // }
  if (isBoundary(ix, iy, iz) ) {
    const SurfaceProperty& sp = surfacePropertyAt(ix, iy, iz);
    const AirProperty& ap1 = innerPropertyAt(ix, iy, iz);
    double twall = sp.T();
    value = ap1.T();
    // std::cout << "Boundary at (" << ix << ", " << iy << ", " << iz
    // 	      << ") twall=" << twall 
    // 	      << " BC=" << sp.boundaryCondition().flags() << std::endl;
    if (sp.hasBoundaryCondition(BoundaryCondition::kConstantTheta) ) {
      // std::cout << "Constant T wall at (" << ix << ", " << iy << ", " << iz
      // 		<< ")" << std::endl;
      value = 2*twall - ap1.T();
    }
  } else {
    const AirProperty& ap0 = mAtmosphere.propertyAt(ix, iy, iz);
    value = ap0.T();
  }
  // ix = (ix + mNPoints[0])%mNPoints[0];
  // iy = (iy + mNPoints[1])%mNPoints[1];
  // if (iz >= 0 && iz < mNPoints[2]) {
  //   const AirProperty& p1 = mAtmosphere.propertyAt(ix, iy, iz);
  //   value = p1.T();
  // } else if (iz==-1) {
  //   const AirProperty& p1 = mAtmosphere.propertyAt(ix, iy, 0);
  //   const SurfaceProperty& p2 = mSurfaceState.propertyAt(ix, iy);
  //   value = 2*p2.T()-p1.T();
  // } else if (iz == mNPoints[2]) {
  //   const AirProperty& p1 = mAtmosphere.propertyAt(ix, iy, mNPoints[2]-1);
  //   value = p1.T();
  // }
  return value;
}

double Environment::u(int ix, int iy, int iz) const {
  double value=0.0;
  if (ix < -1) {
    value = 0.0;
  } else if (isBoundary(ix, iy, iz) ) {
    int naxis = 0;
    if (ix < 0 || ix >= mNPoints[0]) naxis ++;
    if (iy < 0 || iy >= mNPoints[1]) naxis ++;
    if (iz < 0 || iz >= mNPoints[2]) naxis ++;
    if (naxis > 1) {
      value = 0.0;
    } else {
      // BoundaryCell bcell = boundaryCell(this, ix, iy, iz);
      // value = bcell.u();
      const SurfaceProperty& sp = surfacePropertyAt(ix, iy, iz);
      const AirProperty& ap1 = innerPropertyAt(ix, iy, iz);
      BoundaryCell::Location loc = surfaceLocation(ix, iy, iz);
      if (loc==BoundaryCell::kXLow || loc==BoundaryCell::kXHigh) {
	value = 0.0;
      } else {
	value = -(ap1.u() );
      }
    }
  } else {
    const AirProperty& ap0 = mAtmosphere.propertyAt(ix, iy, iz);
    value = ap0.u();
    if (ix == (mNPoints[0]-1) ) {
      value = 0.0;
    }
  }
  // ix = (ix + mNPoints[0])%mNPoints[0];
  // iy = (iy + mNPoints[1])%mNPoints[1];
  // if (iz >= 0 && iz < mNPoints[2]) {
  //   const AirProperty& p1 = mAtmosphere.propertyAt(ix, iy, iz);
  //   value = p1.u();
  // } else if (iz==-1) {
  //   const AirProperty& p1 = mAtmosphere.propertyAt(ix, iy, 0);
  //   value = -p1.u();
  // } else if (iz == mNPoints[2]) {
  //   const AirProperty& p1 = mAtmosphere.propertyAt(ix, iy, mNPoints[2]-1);
  //   value = -p1.u();
  // }
  return value;
}

double Environment::v(int ix, int iy, int iz) const {
  double value=0.0;
  if (iy < -1) {
    value = 0.0;
  } else if (isBoundary(ix, iy, iz) ) {
    int naxis = 0;
    if (ix < 0 || ix >= mNPoints[0]) naxis ++;
    if (iy < 0 || iy >= mNPoints[1]) naxis ++;
    if (iz < 0 || iz >= mNPoints[2]) naxis ++;
    if (naxis > 1) {
      value = 0.0;
    } else {
      const SurfaceProperty& sp = surfacePropertyAt(ix, iy, iz);
      const AirProperty& ap1 = innerPropertyAt(ix, iy, iz);
      BoundaryCell::Location loc = surfaceLocation(ix, iy, iz);
      if (loc==BoundaryCell::kYLow || loc==BoundaryCell::kYHigh) {
	value = 0.0;
      } else {
	value = -(ap1.v() );
      }
    }
  } else {
    const AirProperty& ap0 = mAtmosphere.propertyAt(ix, iy, iz);
    value = ap0.v();
    if (iy == (mNPoints[1]-1) ) {
      value = 0.0;
    }
  }
  // ix = (ix + mNPoints[0])%mNPoints[0];
  // iy = (iy + mNPoints[1])%mNPoints[1];
  // if (iz >= 0 && iz < mNPoints[2]) {
  //   const AirProperty& p1 = mAtmosphere.propertyAt(ix, iy, iz);
  //   value = p1.v();
  // } else if (iz==-1) {
  //   const AirProperty& p1 = mAtmosphere.propertyAt(ix, iy, 0);
  //   value = -p1.v();
  // } else if (iz == mNPoints[2]) {
  //   const AirProperty& p1 = mAtmosphere.propertyAt(ix, iy, mNPoints[2]-1);
  //   value = -p1.v();
  // }
  return value;
}

double Environment::w(int ix, int iy, int iz) const {
  double value=0.0;
  if (iz < -1) {
  } else if (isBoundary(ix, iy, iz) ) {
    int naxis = 0;
    if (ix < 0 || ix >= mNPoints[0]) naxis ++;
    if (iy < 0 || iy >= mNPoints[1]) naxis ++;
    if (iz < 0 || iz >= mNPoints[2]) naxis ++;
    if (naxis > 1) {
      value = 0.0;
    } else {
      const SurfaceProperty& sp = surfacePropertyAt(ix, iy, iz);
      const AirProperty& ap1 = innerPropertyAt(ix, iy, iz);
      BoundaryCell::Location loc = surfaceLocation(ix, iy, iz);
      if (loc==BoundaryCell::kZLow || loc==BoundaryCell::kZHigh) {
	value = 0.0;
      } else {
	value = -(ap1.w() );
      }
    }
  } else {
    const AirProperty& ap0 = mAtmosphere.propertyAt(ix, iy, iz);
    value = ap0.w();
    if (iy == (mNPoints[2]-1) ) {
      value = 0.0;
    }
  }
  // ix = (ix + mNPoints[0])%mNPoints[0];
  // iy = (iy + mNPoints[1])%mNPoints[1];
  // if (iz >= 0 && iz < (mNPoints[2]-1) ) {
  //   const AirProperty& p1 = mAtmosphere.propertyAt(ix, iy, iz);
  //   value = p1.w();
  // } else if (iz==-1) {
  //   value = 0.0;
  // } else if (iz >= (mNPoints[2]-1) ) {
  //   value = 0.0;
  // }
  return value;
}

bool Environment::isBoundary(int ix, int iy, int iz) const {
  bool r=false;
  bool xin = (ix >= 0 && ix < mNPoints[0]);
  bool yin = (iy >= 0 && iy < mNPoints[1]);
  bool zin = (iz >= 0 && iz < mNPoints[2]);
  return !(xin && yin && zin);
}

BoundaryCell::Location 
Environment::surfaceLocation(int ix, int iy, int iz) const {
  BoundaryCell::Location loc;
  if (ix < 0) {
    loc = BoundaryCell::kXLow;
  } else if (ix >= mNPoints[0]) {
    loc = BoundaryCell::kXHigh;
  } else if (iy < 0) {
    loc = BoundaryCell::kYLow;
  } else if (iy >= mNPoints[1]) {
    loc = BoundaryCell::kYHigh;
  } else if (iz < 0) {
    loc = BoundaryCell::kZLow;
  } else if (iz >= mNPoints[2]) {
    loc = BoundaryCell::kZHigh;
  }
  return loc;
}

const SurfaceProperty& 
Environment::surfacePropertyAt(int ix, int iy, int iz) const {
  SurfaceProperty spsp;
  const SurfaceProperty& sp=spsp;

   // std::cout << "(ix, iy, iz)=(" << ix << ", " << iy << ", " << iz << ")"
   // 	    << std::endl;
  if (ix < 0) {
    //    std::cout << "ix<0 matched" << std::endl;
    if (mSurfaceStates[BoundaryCell::kXLow]) {
      //      std::cout << "return surfaceState XLow" << std::endl;
      //      mSurfaceStates[BoundaryCell::kXLow]->dump();
      //      mSurfaceStates[BoundaryCell::kXLow]->propertyAt(iy, iz).print();
      return mSurfaceStates[BoundaryCell::kXLow]->propertyAt(iy, iz);
    }
    //    std::cout << "ix<0 through" << std::endl;
  } else if (ix >= mNPoints[0]) {
    //    std::cout << "ix>=nx matched" << std::endl;
    if (mSurfaceStates[BoundaryCell::kXHigh]) {
      //      std::cout << "return surfaceState XHigh" << std::endl;
      return mSurfaceStates[BoundaryCell::kXHigh]->propertyAt(iy, iz);
    }
  } else if (iy < 0) {
    if (mSurfaceStates[BoundaryCell::kYLow]) {
      return mSurfaceStates[BoundaryCell::kYLow]->propertyAt(iz, ix);
    }
  } else if (iy >= mNPoints[1]) {
    if (mSurfaceStates[BoundaryCell::kYHigh]) {
      return mSurfaceStates[BoundaryCell::kYHigh]->propertyAt(iz, ix);
    }
  } else if (iz < 0) {
    if (mSurfaceStates[BoundaryCell::kZLow]) {
      return mSurfaceStates[BoundaryCell::kZLow]->propertyAt(ix, iy);
    }
  } else if (iz >= mNPoints[2]) {
    if (mSurfaceStates[BoundaryCell::kZHigh]) {
      return mSurfaceStates[BoundaryCell::kZHigh]->propertyAt(ix, iy);
    }
  }
  return sp;
}

const AirProperty& Environment::innerPropertyAt(int ix, int iy, int iz) const {
  // std::cout << "Get inner property at (" << ix << ", " << iy << ", " << iz 
  // 	    << ")" << std::endl;
  if (ix < 0) {
    return mAtmosphere.propertyAt(ix+1, iy, iz);
  } else if (ix >= mNPoints[0]) {
    return mAtmosphere.propertyAt(ix-1, iy, iz);
  } if (iy < 0) {
    return mAtmosphere.propertyAt(ix, iy+1, iz);
  } else if (iy >= mNPoints[1]) {
    return mAtmosphere.propertyAt(ix, iy-1, iz);
  } if (iz < 0) {
    return mAtmosphere.propertyAt(ix, iy, iz+1);
  } else if (iz >= mNPoints[2]) {
    return mAtmosphere.propertyAt(ix, iy, iz-1);
  }
  return mAtmosphere.propertyAt(ix, iy, iz);
}

void Environment::dump() {
}


