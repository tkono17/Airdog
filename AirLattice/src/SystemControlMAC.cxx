/*
  SystemControlMAC.cxx
*/
#include <iostream>
#include "AirLattice/SystemControl.hxx"
#include "AirLattice/SystemControlMAC.hxx"
#include "AirLattice/ScaleSet.hxx"
#include "AirLattice/LatticeIndexTool.hxx"
#include "AirLattice/MatrixTool.hxx"
#include "AirLattice/Atmosphere.hxx"
#include "AirLattice/AirProperty.hxx"
#include "AirLattice/Constants.hxx"

SystemControlMAC::SystemControlMAC() {
}

SystemControlMAC::SystemControlMAC(Environment* env, double dt, double ntimepoints) {
}

SystemControlMAC::~SystemControlMAC() {
}

int SystemControlMAC::initialize() {
  SystemControl::initialize();
  mLatticeIndexTool.initialize(mEnvironment);
  uniformHotspot2();

  ScaleSet::get()->print();

  return 0;
}

int SystemControlMAC::run() {
  return SystemControl::run();
}

void SystemControlMAC::updateEnvironment(float dt) {
  Atmosphere& atmosphere = mEnvironment->atmosphere();
  // std::cout << "Field at (ix, iy, iz)=(49-51, 0, 1)" << std::endl;
  // atmosphere.propertyAt(49, 0, 1).print();
  // atmosphere.propertyAt(50, 0, 1).print();
  // atmosphere.propertyAt(51, 0, 1).print();
  int i1=0, i2=1, i3=2;
  // std::cout << "Field at (ix, iy, iz)=(49-51, 0, 0)" << std::endl;
  // atmosphere.propertyAt(i1, 0, 0).print();
  // atmosphere.propertyAt(i2, 0, 0).print();
  // atmosphere.propertyAt(i3, 0, 0).print();

  int ncells = mLatticeIndexTool.NCells();
  mM.clear();
  mM.initialize(ncells, ncells);
  mB.assign(ncells, 0.0);

  //  std::cout << "Update pressure (ncells=" << ncells << ")" << std::endl;
  loop(std::mem_fun(&SystemControlMAC::buildMatrixForP) );
  // printMatrix(mM);
  // printVector(mB);
  //  std::cout << "Solve equation to obtain P" << std::endl;
  solveP();

  //  std::cout << "Update velocity" << std::endl;
  loop(std::mem_fun(&SystemControlMAC::updateVelocity) );
  //  std::cout << "  -> update variables" << std::endl;
  loop(std::mem_fun(&SystemControlMAC::updateVariables) );

  //  std::cout << "Update temperature" << std::endl;
  loop(std::mem_fun(&SystemControlMAC::updateTemperature) );
  //  std::cout << "  -> update variables" << std::endl;
  loop(std::mem_fun(&SystemControlMAC::updateVariables) );

}

void SystemControlMAC::applyInitialConditions() {
  return uniformHotspot2();
}

void SystemControlMAC::updateTemperature() {
  //  SystemControl::updateTemperature();
  int ix = mXYZIndex[0];
  int iy = mXYZIndex[1];
  int iz = mXYZIndex[2];
  //  ScaleSet* ss = ScaleSet::get();

  double t = mDeltaTime*(-convTheta(ix, iy, iz) + laplacianTheta(ix, iy, iz) );

  Atmosphere& atmosphere = mEnvironment->atmosphere();
  AirProperty& ap = atmosphere.propertyAt(ix, iy, iz);
  ap.updateVariable(kT, t);
}

void SystemControlMAC::buildMatrixForP() {
  int ix=mXYZIndex[0];
  int iy=mXYZIndex[1];
  int iz=mXYZIndex[2];

  ScaleSet* ss = ScaleSet::get();

  int ncells = mLatticeIndexTool.NCells();

  // Update matrix M and vector b (Mx = b) around this cell
  // Matrix M (affects M 1+6 elements)
  //------------------------------------------------------------
  int nz = mLatticeIndexTool.Nz();
  double dx2=mDeltaX*mDeltaX;
  double dy2=mDeltaY*mDeltaY;
  double dz2=mDeltaZ*mDeltaZ;
  int index0 = mLatticeIndexTool.indexOfCell(ix, iy, iz);
  int index_x1 = mLatticeIndexTool.indexOfCell(ix-1, iy, iz);
  int index_x2 = mLatticeIndexTool.indexOfCell(ix+1, iy, iz);
  int index_y1 = mLatticeIndexTool.indexOfCell(ix, iy-1, iz);
  int index_y2 = mLatticeIndexTool.indexOfCell(ix, iy+1, iz);
  int index_z1 = mLatticeIndexTool.indexOfCell(ix, iy, iz-1);
  int index_z2 = mLatticeIndexTool.indexOfCell(ix, iy, iz+1);
  double value0=0.0;
  // x1, x2
  mM.setValue(index0, index_x1, 1/dx2);
  mM.setValue(index0, index_x2, 1/dx2);
  value0 += -2/dx2;
  // y1, y2
  mM.setValue(index0, index_y1, 1/dy2);
  mM.setValue(index0, index_y2, 1/dy2);
  value0 += -2/dy2;
  // z1, z2 and this cell (ijk)
  if (iz == 0) {
    value0 += -1/dz2;
    mM.setValue(index0, index_z2, 1/dz2);
    mM.setValue(index0, index0, value0);
  } else if (iz == (nz-1) ) {
    value0 += -1/dz2;
    mM.setValue(index0, index_z1, 1/dz2);
    mM.setValue(index0, index0, value0);
  } else {
    value0 += -2/dz2;
    mM.setValue(index0, index_z1, 1/dz2);
    mM.setValue(index0, index_z2, 1/dz2);
    mM.setValue(index0, index0, value0);
  }

  // Vector b (affects this index only)
  //------------------------------------------------------------
  double b=0.0;
  double pr = ss->Pr();
  double bc = ss->buoyancyCoefficient();

  //  std::cout << " pr=" << pr << ", bc=" << bc << std::endl;
 
  b += divV(ix, iy, iz)/mDeltaTime;
  b += ( (pr*laplacianU(ix, iy, iz) - convU(ix, iy, iz) )
	 - (pr*laplacianU(ix-1, iy, iz) - convU(ix-1, iy, iz) ) )/mDeltaX;
  b += ( (pr*laplacianV(ix, iy, iz) - convV(ix, iy, iz) )
	 - (pr*laplacianV(ix, iy-1, iz) - convV(ix, iy-1, iz) ) )/mDeltaY;
  b += ( (pr*laplacianW(ix, iy, iz) - convW(ix, iy, iz) )
	 - (pr*laplacianW(ix, iy, iz-1) - convW(ix, iy, iz-1) ) )/mDeltaZ;
  b += (buoyancy(ix, iy, iz)-buoyancy(ix, iy, iz-1) )/mDeltaZ;
  mB[index0] = b;
}

void SystemControlMAC::solveP() {
  solveGausSeidel(mM, mB);
  unsigned int n=mB.size();
  unsigned int i;
  Atmosphere& atmosphere = mEnvironment->atmosphere();
  std::vector<int> v;

  for (i=0; i<n; ++i) {
    v = mLatticeIndexTool.cellIndexAt(i);
    AirProperty& ap = atmosphere.propertyAt(v[0], v[1], v[2]);
    ap.setPressure(mB[i]);
  }
}

void SystemControlMAC::updatePressure() {
  SystemControl::updatePressure();
}

void SystemControlMAC::updateVelocity() {
  //  SystemControl::updateVelocity();
  double v[3] = { 0, 0, 0 };
  int ix = mXYZIndex[0];
  int iy = mXYZIndex[1];
  int iz = mXYZIndex[2];
  ScaleSet* ss = ScaleSet::get();
  double pr = ss->Pr();

  double dpdx=(mEnvironment->P(ix+1, iy, iz) - mEnvironment->P(ix, iy, iz) )/mDeltaX;
  double dpdy=(mEnvironment->P(ix, iy+1, iz) - mEnvironment->P(ix, iy, iz) )/mDeltaY;
  double dpdz=(mEnvironment->P(ix, iy, iz+1) - mEnvironment->P(ix, iy, iz) )/mDeltaZ;
  v[0] = mDeltaTime*(-convU(ix, iy, iz) - dpdx + pr*laplacianU(ix, iy, iz) );
  v[1] = mDeltaTime*(-convV(ix, iy, iz) - dpdy + pr*laplacianV(ix, iy, iz) );
  v[2] = mDeltaTime*(-convW(ix, iy, iz) - dpdz + pr*laplacianW(ix, iy, iz) );

  Atmosphere& atmosphere = mEnvironment->atmosphere();
  AirProperty& ap = atmosphere.propertyAt(ix, iy, iz);
  ap.updateVariable(kU, v[0]);
  ap.updateVariable(kV, v[1]);
  ap.updateVariable(kW, v[2]);
}

void SystemControlMAC::applyBoundaryConditions() {
  SystemControl::applyBoundaryConditions();
}

void SystemControlMAC::updateVariables() {
  SystemControl::updateVariables();
}

double SystemControlMAC::convU(int ix, int iy, int iz) {
  double y=0.0;
  double u=mEnvironment->u(ix, iy, iz);
  double v=(mEnvironment->v(ix,iy,iz) + mEnvironment->v(ix,iy-1,iz) +
	    mEnvironment->v(ix+1,iy,iz) + mEnvironment->v(ix+1,iy-1,iz) )/4;
  double w=(mEnvironment->w(ix,iy,iz) + mEnvironment->w(ix,iy,iz-1) +
	    mEnvironment->w(ix+1,iy,iz) + mEnvironment->w(ix+1,iy,iz-1) )/4;

  double u_0, u_m1, u_p1;
  u_0 = mEnvironment->u(ix, iy, iz);

  u_m1 = mEnvironment->u(ix-1, iy, iz);
  u_p1 = mEnvironment->u(ix+1, iy, iz);
  y += u*(u_p1 - u_m1)/(2*mDeltaX) 
    - std::fabs(u)*(u_p1 - 2*u_0 + u_m1)/(2*mDeltaX);

  u_m1 = mEnvironment->u(ix, iy-1, iz);
  u_p1 = mEnvironment->u(ix, iy+1, iz);
  y += v*(u_p1 - u_m1)/(2*mDeltaY) 
    - std::fabs(v)*(u_p1 - 2*u_0 + u_m1)/(2*mDeltaY);

  u_m1 = mEnvironment->u(ix, iy, iz-1);
  u_p1 = mEnvironment->u(ix, iy, iz+1);
  y += w*(u_p1 - u_m1)/(2*mDeltaZ) 
    - std::fabs(w)*(u_p1 - 2*u_0 + u_m1)/(2*mDeltaZ);
  
  return y;
}

double SystemControlMAC::convV(int ix, int iy, int iz) {
  double y=0.0;
  double u=(mEnvironment->u(ix,iy,iz) + mEnvironment->v(ix-1,iy,iz) +
	    mEnvironment->u(ix,iy+1,iz) + mEnvironment->v(ix-1,iy+1,iz) )/4;
  double v=mEnvironment->v(ix, iy, iz);
  double w=(mEnvironment->w(ix,iy,iz) + mEnvironment->w(ix,iy,iz-1) +
	    mEnvironment->w(ix,iy+1,iz) + mEnvironment->w(ix,iy+1,iz-1) )/4;

  double v_0, v_m1, v_p1;
  v_0 = mEnvironment->v(ix, iy, iz);

  v_m1 = mEnvironment->v(ix-1, iy, iz);
  v_p1 = mEnvironment->v(ix+1, iy, iz);
  y += u*(v_p1 - v_m1)/(2*mDeltaX) 
    - std::fabs(u)*(v_p1 - 2*v_0 + v_m1)/(2*mDeltaX);

  v_m1 = mEnvironment->v(ix, iy-1, iz);
  v_p1 = mEnvironment->v(ix, iy+1, iz);
  y += v*(v_p1 - v_m1)/(2*mDeltaY) 
    - std::fabs(v)*(v_p1 - 2*v_0 + v_m1)/(2*mDeltaY);

  v_m1 = mEnvironment->v(ix, iy, iz-1);
  v_p1 = mEnvironment->v(ix, iy, iz+1);
  y += w*(v_p1 - v_m1)/(2*mDeltaZ) 
    - std::fabs(w)*(v_p1 - 2*v_0 + v_m1)/(2*mDeltaZ);
  
  return y;
}

double SystemControlMAC::convW(int ix, int iy, int iz) {
  double y=0.0;
  double u=mEnvironment->u(ix, iy, iz);
  double v=(mEnvironment->v(ix,iy,iz) + mEnvironment->v(ix,iy-1,iz) +
	    mEnvironment->v(ix+1,iy,iz) + mEnvironment->v(ix+1,iy-1,iz) )/4;
  double w=(mEnvironment->w(ix,iy,iz) + mEnvironment->w(ix,iy,iz-1) +
	    mEnvironment->w(ix+1,iy,iz) + mEnvironment->w(ix+1,iy,iz-1) )/4;

  double w_0, w_m1, w_p1;
  w_0 = mEnvironment->w(ix+1, iy, iz);

  w_m1 = mEnvironment->w(ix-1, iy, iz);
  w_p1 = mEnvironment->w(ix+1, iy, iz);
  y += u*(w_p1 - w_m1)/(2*mDeltaX) 
    - std::fabs(u)*(w_p1 - 2*w_0 + w_m1)/(2*mDeltaX);

  w_m1 = mEnvironment->w(ix, iy-1, iz);
  w_p1 = mEnvironment->w(ix, iy+1, iz);
  y += v*(w_p1 - w_m1)/(2*mDeltaY) 
    - std::fabs(v)*(w_p1 - 2*w_0 + w_m1)/(2*mDeltaY);

  w_m1 = mEnvironment->w(ix, iy, iz-1);
  w_p1 = mEnvironment->w(ix, iy, iz+1);
  y += w*(w_p1 - w_m1)/(2*mDeltaZ) 
    - std::fabs(w)*(w_p1 - 2*w_0 + w_m1)/(2*mDeltaZ);
  
  return y;
}

double SystemControlMAC::convTheta(int ix, int iy, int iz) {
  double y=0.0;
  double u=(mEnvironment->u(ix, iy, iz) + mEnvironment->u(ix-1, iy, iz) )/2;
  double v=(mEnvironment->v(ix, iy, iz) + mEnvironment->v(ix, iy-1, iz) )/2;
  double w=(mEnvironment->w(ix, iy, iz) + mEnvironment->w(ix, iy, iz-1) )/2;

  double t_0, t_m1, t_p1;
  t_0 = mEnvironment->theta(ix, iy, iz);

  t_m1 = mEnvironment->theta(ix-1, iy, iz);
  t_p1 = mEnvironment->theta(ix+1, iy, iz);
  y += u*(t_p1 - t_m1)/(2*mDeltaX) 
    - std::fabs(u)*(t_p1 - 2*t_0 + t_m1)/(2*mDeltaX);

  t_m1 = mEnvironment->theta(ix, iy-1, iz);
  t_p1 = mEnvironment->theta(ix, iy+1, iz);
  y += v*(t_p1 - t_m1)/(2*mDeltaY) 
    - std::fabs(v)*(t_p1 - 2*t_0 + t_m1)/(2*mDeltaY);

  t_m1 = mEnvironment->theta(ix, iy, iz-1);
  t_p1 = mEnvironment->theta(ix, iy, iz+1);
  y += w*(t_p1 - t_m1)/(2*mDeltaZ) 
    - std::fabs(w)*(t_p1 - 2*t_0 + t_m1)/(2*mDeltaZ);
  
  return y;
}

double SystemControlMAC::divV(int ix, int iy, int iz) {
  double y=0.0;
  y += (mEnvironment->u(ix,iy,iz)-mEnvironment->u(ix-1,iy,iz) )/mDeltaX;
  y += (mEnvironment->v(ix,iy,iz)-mEnvironment->v(ix,iy-1,iz) )/mDeltaY;
  y += (mEnvironment->w(ix,iy,iz)-mEnvironment->w(ix,iy,iz-1) )/mDeltaZ;
  return y;
}

double SystemControlMAC::buoyancy(int ix, int iy, int iz) {
  double c=1.0;
  ScaleSet* ss = ScaleSet::get();
  if (ss) {
    c = ss->buoyancyCoefficient();
  }
  double t=0.0;
  t = (mEnvironment->theta(ix, iy, iz) + mEnvironment->theta(ix, iy, iz+1) )/2;
  return c*t;
}

double SystemControlMAC::laplacianU(int ix, int iy, int iz) {
  double f0=mEnvironment->u(ix, iy, iz);
  double fx1=mEnvironment->u(ix-1, iy, iz);
  double fx2=mEnvironment->u(ix+1, iy, iz);
  double fy1=mEnvironment->u(ix, iy-1, iz);
  double fy2=mEnvironment->u(ix, iy+1, iz);
  double fz1=mEnvironment->u(ix, iy, iz-1);
  double fz2=mEnvironment->u(ix, iy, iz+1);
  double y=0.0;
  y += (fx2 - 2*f0 + fx1)/(mDeltaX*mDeltaX) 
    + (fy2 - 2*f0 + fy1)/(mDeltaY*mDeltaY) 
    + (fz2 - 2*f0 + fz1)/(mDeltaZ*mDeltaZ);
  return y;
}

double SystemControlMAC::laplacianV(int ix, int iy, int iz) {
  double f0=mEnvironment->v(ix, iy, iz);
  double fx1=mEnvironment->v(ix-1, iy, iz);
  double fx2=mEnvironment->v(ix+1, iy, iz);
  double fy1=mEnvironment->v(ix, iy-1, iz);
  double fy2=mEnvironment->v(ix, iy+1, iz);
  double fz1=mEnvironment->v(ix, iy, iz-1);
  double fz2=mEnvironment->v(ix, iy, iz+1);
  double y=0.0;
  y += (fx2 - 2*f0 + fx1)/(mDeltaX*mDeltaX) 
    + (fy2 - 2*f0 + fy1)/(mDeltaY*mDeltaY) 
    + (fz2 - 2*f0 + fz1)/(mDeltaZ*mDeltaZ);
  return y;
}

double SystemControlMAC::laplacianW(int ix, int iy, int iz) {
  double f0=mEnvironment->w(ix, iy, iz);
  double fx1=mEnvironment->w(ix-1, iy, iz);
  double fx2=mEnvironment->w(ix+1, iy, iz);
  double fy1=mEnvironment->w(ix, iy-1, iz);
  double fy2=mEnvironment->w(ix, iy+1, iz);
  double fz1=mEnvironment->w(ix, iy, iz-1);
  double fz2=mEnvironment->w(ix, iy, iz+1);
  double y=0.0;
  y += (fx2 - 2*f0 + fx1)/(mDeltaX*mDeltaX) 
    + (fy2 - 2*f0 + fy1)/(mDeltaY*mDeltaY) 
    + (fz2 - 2*f0 + fz1)/(mDeltaZ*mDeltaZ);
  return y;
}

double SystemControlMAC::laplacianTheta(int ix, int iy, int iz) {
  double f0=mEnvironment->theta(ix, iy, iz);
  double fx1=mEnvironment->theta(ix-1, iy, iz);
  double fx2=mEnvironment->theta(ix+1, iy, iz);
  double fy1=mEnvironment->theta(ix, iy-1, iz);
  double fy2=mEnvironment->theta(ix, iy+1, iz);
  double fz1=mEnvironment->theta(ix, iy, iz-1);
  double fz2=mEnvironment->theta(ix, iy, iz+1);
  double y=0.0;
  y += (fx2-2*f0+fx1)/(mDeltaX*mDeltaX) 
    + (fy2-2*f0+fy1)/(mDeltaY*mDeltaY) 
    + (fz2-2*f0+fz1)/(mDeltaZ*mDeltaZ);
  return y;
}

void SystemControlMAC::loop(std::mem_fun_t<void, SystemControlMAC> action) {
  int ix, iy, iz;
  int ix1, ix2, iy1, iy2, iz1, iz2;
  const int* ns = mEnvironment->NPoints();

  Atmosphere& atmosphere = mEnvironment->atmosphere();
  //  SurfaceState& surface = mEnvironment->surfaceState();

  for (ix=0; ix<ns[0]; ++ix) {
    for (iy=0; iy<ns[1]; ++iy) {
      for (iz=0; iz<ns[2]; ++iz) {
	ix1 = ix - 1;
	ix2 = ix + 1;
	iy1 = iy - 1;
	iy2 = iy + 1;
	iz1 = iz - 1;
	iz2 = iz + 1;

	mXYZIndex[0] = ix;
	mXYZIndex[1] = iy;
	mXYZIndex[2] = iz;

	if ( (mXYZIndex[0]==50 || mXYZIndex[0]==49 || mXYZIndex[0]==51)
	    && mXYZIndex[1]==0 && mXYZIndex[2]==0) {
	  mDoPrint = true;
	} else {
	  mDoPrint = false;
	}
	mAP0 = &atmosphere.propertyAt(ix, iy, iz);
	// mAPx1 = &atmosphere.propertyAt(ix1, iy, iz);
	// mAPx2 = &atmosphere.propertyAt(ix2, iy, iz);
	// mAPy1 = &atmosphere.propertyAt(ix, iy1, iz);
	// mAPy2 = &atmosphere.propertyAt(ix, iy2, iz);
	// if (iz1 < 0) {
	//   mAPz1 = 0;
	// } else {
	//   mAPz1 = &atmosphere.propertyAt(ix, iy, iz1);
	// }
	// if (iz2 < 0) {
	//   mAPz2 = 0;
	// } else {
	//   mAPz2 = &atmosphere.propertyAt(ix, iy, iz2);
	// }

	action(this);
      }
    }
  }
}

void SystemControlMAC::uniformHotspot2() {
  int ix, iy, iz;
  const int* np = mEnvironment->NPoints();
  Atmosphere& atmosphere = mEnvironment->atmosphere();
  SurfaceState& surface = mEnvironment->surfaceState();
  const double t0 = Constants::T0;

  std::cout << "UniformHotspot2 condition" << std::endl;

  int nx = np[0];
  int ny = np[1];
  int nz = np[2];
  int hsx = static_cast<int>(np[0]/2);
  int hsy = static_cast<int>(np[1]/2);
  double z=0.0;
  double dt = 200.0/nz;

  // std::cout << "  Create a hotspot on the surface ("
  // 	    << hsx << ", " << hsy << ")" << std::endl;

  ScaleSet* ss = ScaleSet::get();
  double rho = ss->rho();

  std::cout << "Lattice : (nx, ny, nz)=(" << nx << ", " << ny << ", " << nz << ")" << std::endl;
  for (ix=0; ix<nx; ++ix) {
    for (iy=0; iy<ny; ++iy) {
      for (iz=0; iz<nz; ++iz) {
	z = mEnvironment->z(iz);
	AirProperty& ap = atmosphere.propertyAt(ix, iy, iz);
	ap.setPressure(1.0);
	ap.setTemperature(0.0);
	ap.setVelocity(0.0, 0.0, 0.0);
	ap.setDensity(rho);
	ap.setVaporPressure(0.0);
	ap.setWaterDensity(0.0);

      }
      SurfaceProperty& sp = surface.propertyAt(ix, iy);
      sp.setTemperature(0);
      if (ix < nx/3) {
	sp.setTemperature(20);
      }
    }
  }
}

