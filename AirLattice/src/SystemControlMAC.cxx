/*
  SystemControlMAC.cxx
*/
#include <iostream>
#include <bitset>
#include "AirLattice/SystemControl.hxx"
#include "AirLattice/SystemControlMAC.hxx"
#include "AirLattice/ScaleSet.hxx"
#include "AirLattice/LatticeIndexTool.hxx"
#include "AirLattice/MatrixTool.hxx"
#include "AirLattice/Atmosphere.hxx"
#include "AirLattice/AirProperty.hxx"
#include "AirLattice/Constants.hxx"
#include "TRandom3.h"

SystemControlMAC::SystemControlMAC() : 
  SystemControl(), mRandom() {
}

SystemControlMAC::SystemControlMAC(Environment* env, double dt, double ntimepoints) : 
  SystemControl(env, dt, ntimepoints), 
  mRandom() {
}

SystemControlMAC::~SystemControlMAC() {
}

int SystemControlMAC::initialize() {
  SystemControl::initialize();
  mLatticeIndexTool.initialize(mEnvironment);
  uniformHotspot2();

  ScaleSet::get()->print();
  std::cout << "(dT, dX, dY, dZ) = (" << mDeltaTime
	    << ", " << mDeltaX
	    << ", " << mDeltaY
	    << ", " << mDeltaZ << ")" << std::endl;
  std::cout << "(dT/dX2, dT/dY2, dT/dZ2) = (" 
	    << (mDeltaTime/(mDeltaX*mDeltaX) ) << ", "
	    << (mDeltaTime/(mDeltaY*mDeltaY) ) << ", "
	    << (mDeltaTime/(mDeltaZ*mDeltaZ) ) << ")"
	    << std::endl;

  mTimeStep = 0;

  int ncells = mLatticeIndexTool.NCells();
  mLinearSolver.assignMatrix(ncells, ncells);
  mLinearSolver.assignVector(ncells);
  mM.clear();
  mM.initialize(ncells, ncells);
  mB.assign(ncells, 0.0);
  mPVec.assign(ncells, 1.0);

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

  if (mTimeStep == 100 && mRecordInterval < 50) {
    mRecordInterval = 50;
  } else if (mTimeStep == 1000 && mRecordInterval < 500) {
    mRecordInterval = 500;
  } else if (mTimeStep == 10000 && mRecordInterval < 5000) {
    mRecordInterval = 5000;
  } else if (mTimeStep == 100000 && mRecordInterval < 50000) {
    mRecordInterval = 50000;
  }

  int ncells = mLatticeIndexTool.NCells();
  for (int i=0; i<ncells; ++i) {
    mLinearSolver.setVectorElement(i, 0.0);
  }
  mB.assign(ncells, 0.0);
  //  std::cout << "Update pressure (ncells=" << ncells << ")" << std::endl;
  loop(std::mem_fun(&SystemControlMAC::buildMatrixForP) );
  //  printMatrix(mM);
  //  printVector(mB);
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

  mTimeStep ++;
}

void SystemControlMAC::applyInitialConditions() {
  return uniformHotspot2();
}

void SystemControlMAC::updateTemperature() {
  //  SystemControl::updateTemperature();
  int ix = mXYZIndex[0];
  int iy = mXYZIndex[1];
  int iz = mXYZIndex[2];
  ScaleSet* ss = ScaleSet::get();
  double a = ss->parameterA();

  a = 1.0; // Pr*Re = 1
  double dt = mDeltaTime*(-convTheta(ix, iy, iz) 
			 + a*laplacianTheta(ix, iy, iz) );

  Atmosphere& atmosphere = mEnvironment->atmosphere();
  AirProperty& ap = atmosphere.propertyAt(ix, iy, iz);
  double t0 = ap.T();

  // if (mTimeStep == 0) {
  //   double r=0.0;
  //   r = mRandom.Gaus(0, 0.2);
  //   t += r;
  // }

  ap.updateVariable(kT, t0 + dt);
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
  int nx = mLatticeIndexTool.Nx();
  int ny = mLatticeIndexTool.Ny();
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
  double b=0.0;

  // if (!mEnvironment->isBoundary(ix, iy, iz) ) {
  //   // x1, x2
  //   mM.setValue(index0, index_x1, 1/dx2);
  //   mM.setValue(index0, index_x2, 1/dx2);
  //   value0 += -2/dx2;
  //   // y1, y2
  //   mM.setValue(index0, index_y1, 1/dy2);
  //   mM.setValue(index0, index_y2, 1/dy2);
  //   value0 += -2/dy2;
  //   // z1, z2 and this cell (ijk)
  //   mM.setValue(index0, index_z1, 1/dz2);
  //   mM.setValue(index0, index_z2, 1/dz2);
  //   value0 += -2/dz2;
  //   // This cell
  //   mM.setValue(index0, index0, value0);
  // } else {

  if (mTimeStep == 0) {
    // X boundary
    if (ix == 0 && nx == 1 ) {
    } else if (ix == 0) {
      value0 += -1/dx2;
      mM.setValue(index0, index_x2, 1/dx2);
      mLinearSolver.setMatrixElement(index0, index_x2, 1/dx2);
    } else if (ix == (nx-1) ) {
      value0 += -1/dx2;
      mM.setValue(index0, index_x1, 1/dx2);
      mLinearSolver.setMatrixElement(index0, index_x1, 1/dx2);
    } else {
      value0 += -2/dx2;
      mM.setValue(index0, index_x1, 1/dx2);
      mM.setValue(index0, index_x2, 1/dx2);
      mLinearSolver.setMatrixElement(index0, index_x1, 1/dx2);
      mLinearSolver.setMatrixElement(index0, index_x2, 1/dx2);
    }
    // Y boundary
    if (iy == 0 && ny == 1 ) {
    } else if (iy == 0) {
      value0 += -1/dy2;
      mM.setValue(index0, index_y2, 1/dy2);
      mLinearSolver.setMatrixElement(index0, index_y2, 1/dy2);
    } else if (iy == (ny-1) ) {
      value0 += -1/dy2;
      mM.setValue(index0, index_y1, 1/dy2);
      mLinearSolver.setMatrixElement(index0, index_y1, 1/dy2);
    } else {
      value0 += -2/dy2;
      mM.setValue(index0, index_y1, 1/dy2);
      mM.setValue(index0, index_y2, 1/dy2);
      mLinearSolver.setMatrixElement(index0, index_y1, 1/dy2);
      mLinearSolver.setMatrixElement(index0, index_y2, 1/dy2);
    }
    // Z boundary
    if (iz == 0 && nz == 1 ) {
    } else if (iz == 0) {
      value0 += -1/dz2;
      mM.setValue(index0, index_z2, 1/dz2);
      mLinearSolver.setMatrixElement(index0, index_z2, 1/dz2);
    } else if (iz == (nz-1) ) {
      value0 += -1/dz2;
      mM.setValue(index0, index_z1, 1/dz2);
      mLinearSolver.setMatrixElement(index0, index_z1, 1/dz2);
    } else {
      value0 += -2/dz2;
      mM.setValue(index0, index_z1, 1/dz2);
      mM.setValue(index0, index_z2, 1/dz2);
      mLinearSolver.setMatrixElement(index0, index_z1, 1/dz2);
      mLinearSolver.setMatrixElement(index0, index_z2, 1/dz2);
    }
    mM.setValue(index0, index0, value0);
    mLinearSolver.setMatrixElement(index0, index0, value0);
  }

  // Vector b (affects this index only)
  //------------------------------------------------------------
  double pr = ss->Pr();
  //  double bc = ss->buoyancyCoefficient();

  //  std::cout << " pr=" << pr << ", bc=" << bc << std::endl;
 
  //  std::cout << "b{0} = " << b << std::endl;
  b += divV(ix, iy, iz)/mDeltaTime; // remove for numerical stability
  //  std::cout << "b{1} = " << b << std::endl;
  b += ( (pr*laplacianU(ix, iy, iz) - convU(ix, iy, iz) )
	 - (pr*laplacianU(ix-1, iy, iz) - convU(ix-1, iy, iz) ) )/mDeltaX;
  //  std::cout << "b{2} = " << b << std::endl;
  b += ( (pr*laplacianV(ix, iy, iz) - convV(ix, iy, iz) )
	 - (pr*laplacianV(ix, iy-1, iz) - convV(ix, iy-1, iz) ) )/mDeltaY;
  //  std::cout << "b{3} = " << b << std::endl;
  b += ( (pr*laplacianW(ix, iy, iz) - convW(ix, iy, iz) )
	 - (pr*laplacianW(ix, iy, iz-1) - convW(ix, iy, iz-1) ) )/mDeltaZ;
  //  std::cout << "b{4} = " << b << std::endl;
  b += (buoyancy(ix, iy, iz)-buoyancy(ix, iy, iz-1) )/mDeltaZ;
  //  std::cout << "b{5} = " << b << std::endl;
  mB[index0] = b;
  mLinearSolver.setVectorElement(index0, b);
}

void SystemControlMAC::solveP() {
  //  std::cout << "run Gauss-Seidel method" << std::endl;

  //  solveSOR(mM, mB, mPVec, 1.9);
  solveGausSeidel(mM, mB, mPVec, 1.9);
  //  std::cout << "Done running Gauss-Seidel method" << std::endl;
  unsigned int n=mB.size();
  unsigned int i;
  Atmosphere& atmosphere = mEnvironment->atmosphere();
  std::vector<int> v;
  double b0 = mB[0]-1;

  // mLinearSolver.solve();
  // Eigen::VectorXd& x = mLinearSolver.vectorX();
  // b0 = x(0) - 1.0;

  for (i=0; i<n; ++i) {
    v = mLatticeIndexTool.cellIndexAt(i);
    AirProperty& ap = atmosphere.propertyAt(v[0], v[1], v[2]);
    mB[i] -= b0;
    // x(i) -= b0;
    // mB[i] = x(i);
    ap.setPressure(mB[i]);
    mPVec[i] = 1.0;
  }
}

void SystemControlMAC::updatePressure() {
  SystemControl::updatePressure();
}

void SystemControlMAC::updateVelocity() {
  //  SystemControl::updateVelocity();
  double v[3] = { 0, 0, 0 };
  double dv[3] = { 0, 0, 0 };
  int ix = mXYZIndex[0];
  int iy = mXYZIndex[1];
  int iz = mXYZIndex[2];
  ScaleSet* ss = ScaleSet::get();
  double pr = ss->Pr();
  int nx = mLatticeIndexTool.Nx();
  int ny = mLatticeIndexTool.Ny();
  int nz = mLatticeIndexTool.Nz();
  double dpdx, dpdy, dpdz;

  //  std::cout << "calculate du" << std::endl;
  if (ix == (nx-1) ) {
    dv[0] = 0.0;
  } else {
    dpdx=(mEnvironment->P(ix+1, iy, iz) - mEnvironment->P(ix, iy, iz) )/mDeltaX;
    dv[0] = mDeltaTime*(-convU(ix, iy, iz) - dpdx + pr*laplacianU(ix, iy, iz) );
  }
  //  std::cout << "calculate dv" << std::endl;
  if (iy == (ny - 1) ) {
    dv[1] = 0.0;
  } else {
    dpdy=(mEnvironment->P(ix, iy+1, iz) - mEnvironment->P(ix, iy, iz) )/mDeltaY;
    dv[1] = mDeltaTime*(-convV(ix, iy, iz) - dpdy + pr*laplacianV(ix, iy, iz) );
  }
  //  std::cout << "calculate dw" << std::endl;
  if (iz == (nz-1) ) {
    dv[2] = 0.0;
  } else {
    dpdz=(mEnvironment->P(ix, iy, iz+1) - mEnvironment->P(ix, iy, iz) )/mDeltaZ;
    dv[2] = mDeltaTime*(- convW(ix, iy, iz) 
			- dpdz 
			+ pr*laplacianW(ix, iy, iz)
			+ buoyancy(ix, iy, iz) 
			);
  }
  //  std::cout << "update velocity" << std::endl;

  Atmosphere& atmosphere = mEnvironment->atmosphere();
  AirProperty& ap = atmosphere.propertyAt(ix, iy, iz);
  v[0] = ap.u() + dv[0];
  v[1] = ap.v() + dv[1];
  v[2] = ap.w() + dv[2];
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
  //  std::cout << "convU: top (" << ix << ", " << iy << ", " << iz << ")" << std::endl;
  double y=0.0;
  double u=mEnvironment->u(ix, iy, iz);

  //  std::cout << "convU: (" << ix << ", " << iy << ", " << iz << ")" << std::endl;
  double u_0, u_m1, u_p1;
  const int* np = mEnvironment->NPoints();
  if ( (ix < 0) || ix >= (np[0]-1) ) { // (u, v, w)=0 on the wall
    y = 0.0;
    return y;
  }

  double v=(mEnvironment->v(ix,iy,iz) + mEnvironment->v(ix,iy-1,iz) +
	    mEnvironment->v(ix+1,iy,iz) + mEnvironment->v(ix+1,iy-1,iz) )/4;
  double w=(mEnvironment->w(ix,iy,iz) + mEnvironment->w(ix,iy,iz-1) +
	    mEnvironment->w(ix+1,iy,iz) + mEnvironment->w(ix+1,iy,iz-1) )/4;
  u_0 = mEnvironment->u(ix, iy, iz);

  //  std::cout << "convU x" << std::endl;
  u_m1 = mEnvironment->u(ix-1, iy, iz);
  u_p1 = mEnvironment->u(ix+1, iy, iz);
  y += u*(u_p1 - u_m1)/(2*mDeltaX) 
    - std::fabs(u)*(u_p1 - 2*u_0 + u_m1)/(2*mDeltaX);

  //std::cout << "convU y" << std::endl;
  u_m1 = mEnvironment->u(ix, iy-1, iz);
  u_p1 = mEnvironment->u(ix, iy+1, iz);
  y += v*(u_p1 - u_m1)/(2*mDeltaY) 
    - std::fabs(v)*(u_p1 - 2*u_0 + u_m1)/(2*mDeltaY);

  //std::cout << "convU z" << std::endl;
  u_m1 = mEnvironment->u(ix, iy, iz-1);
  u_p1 = mEnvironment->u(ix, iy, iz+1);
  y += w*(u_p1 - u_m1)/(2*mDeltaZ) 
    - std::fabs(w)*(u_p1 - 2*u_0 + u_m1)/(2*mDeltaZ);
  
  //std::cout << "convU done" << std::endl;
  return y;
}

double SystemControlMAC::convV(int ix, int iy, int iz) {
  double y=0.0;
  double v=mEnvironment->v(ix, iy, iz);

  double v_0, v_m1, v_p1;
  const int* np = mEnvironment->NPoints();
  if ( (iy < 0) || iy >= (np[1]-1) ) {
    y = 0.0;
    return y;
  }

  double u=(mEnvironment->u(ix,iy,iz) + mEnvironment->v(ix-1,iy,iz) +
	    mEnvironment->u(ix,iy+1,iz) + mEnvironment->v(ix-1,iy+1,iz) )/4;
  double w=(mEnvironment->w(ix,iy,iz) + mEnvironment->w(ix,iy,iz-1) +
	    mEnvironment->w(ix,iy+1,iz) + mEnvironment->w(ix,iy+1,iz-1) )/4;
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
  double w=mEnvironment->w(ix, iy, iz);

  double w_0, w_m1, w_p1;
  const int* np = mEnvironment->NPoints();
  if ( (iz < 0) || iz >= (np[2]-1) ) {
    y = 0.0;
    return y;
  }

  double u=(mEnvironment->u(ix,iy,iz) + mEnvironment->u(ix-1,iy,iz) +
	    mEnvironment->u(ix,iy,iz+1) + mEnvironment->u(ix-1,iy,iz+1) )/4;
  double v=(mEnvironment->v(ix,iy,iz) + mEnvironment->v(ix,iy-1,iz) +
	    mEnvironment->v(ix,iy,iz+1) + mEnvironment->v(ix,iy-1,iz+1) )/4;
  w_0 = mEnvironment->w(ix, iy, iz);

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
  //  std::cout << "divV{0}=" << y << std::endl;
  y += (mEnvironment->u(ix,iy,iz)-mEnvironment->u(ix-1,iy,iz) )/mDeltaX;
  //  std::cout << "divV{1}=" << y << std::endl;
  y += (mEnvironment->v(ix,iy,iz)-mEnvironment->v(ix,iy-1,iz) )/mDeltaY;
  //  std::cout << "divV{2}=" << y << std::endl;
  y += (mEnvironment->w(ix,iy,iz)-mEnvironment->w(ix,iy,iz-1) )/mDeltaZ;
  //  std::cout << "divV{3}=" << y << std::endl;
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
  double y=0.0;
  if (ix == -1) {
    double fx2=mEnvironment->u(ix+1, iy, iz);
    y = fx2/(mDeltaX*mDeltaX);
    return y;
  }

  double f0=mEnvironment->u(ix, iy, iz);
  double fx1=mEnvironment->u(ix-1, iy, iz);
  double fx2=mEnvironment->u(ix+1, iy, iz);
  double fy1=mEnvironment->u(ix, iy-1, iz);
  double fy2=mEnvironment->u(ix, iy+1, iz);
  double fz1=mEnvironment->u(ix, iy, iz-1);
  double fz2=mEnvironment->u(ix, iy, iz+1);
  y += (fx2 - 2*f0 + fx1)/(mDeltaX*mDeltaX) 
    + (fy2 - 2*f0 + fy1)/(mDeltaY*mDeltaY) 
    + (fz2 - 2*f0 + fz1)/(mDeltaZ*mDeltaZ);
  return y;
}

double SystemControlMAC::laplacianV(int ix, int iy, int iz) {
  double y=0.0;
  if (iy == -1) {
    double fy2=mEnvironment->v(ix, iy+1, iz);
    y = fy2/(mDeltaY*mDeltaY);
    return y;
  }

  double f0=mEnvironment->v(ix, iy, iz);
  double fx1=mEnvironment->v(ix-1, iy, iz);
  double fx2=mEnvironment->v(ix+1, iy, iz);
  double fy1=mEnvironment->v(ix, iy-1, iz);
  double fy2=mEnvironment->v(ix, iy+1, iz);
  double fz1=mEnvironment->v(ix, iy, iz-1);
  double fz2=mEnvironment->v(ix, iy, iz+1);
  y += (fx2 - 2*f0 + fx1)/(mDeltaX*mDeltaX) 
    + (fy2 - 2*f0 + fy1)/(mDeltaY*mDeltaY) 
    + (fz2 - 2*f0 + fz1)/(mDeltaZ*mDeltaZ);
  return y;
}

double SystemControlMAC::laplacianW(int ix, int iy, int iz) {
  double y=0.0;
  if (iz == -1) {
    double fz2=mEnvironment->w(ix, iy, iz+1);
    y = fz2/(mDeltaZ*mDeltaZ);
    return y;
  }

  double f0=mEnvironment->w(ix, iy, iz);
  double fx1=mEnvironment->w(ix-1, iy, iz);
  double fx2=mEnvironment->w(ix+1, iy, iz);
  double fy1=mEnvironment->w(ix, iy-1, iz);
  double fy2=mEnvironment->w(ix, iy+1, iz);
  double fz1=mEnvironment->w(ix, iy, iz-1);
  double fz2=mEnvironment->w(ix, iy, iz+1);
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
	mXYZIndex[0] = ix;
	mXYZIndex[1] = iy;
	mXYZIndex[2] = iz;

	// std::cout << "Check field property at ("
	// 	  << ix << ", " << iy << ", " << iz << ")"
	// 	  << std::endl;
	if ( (mXYZIndex[0]==50 || mXYZIndex[0]==49 || mXYZIndex[0]==51)
	    && mXYZIndex[1]==0 && mXYZIndex[2]==0) {
	  mDoPrint = true;
	} else {
	  mDoPrint = false;
	}
	mAP0 = &atmosphere.propertyAt(ix, iy, iz);

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
    }
  }

  SurfaceState* surface1=0;
  SurfaceState* surface2=0;

  surface1 = mEnvironment->surfaceState(BoundaryCell::kXLow);
  surface2 = mEnvironment->surfaceState(BoundaryCell::kXHigh);
  for (iy=0; iy<ny; iy++) {
    for (iz=0; iz<nz; ++iz) {
      if (surface1) {
	SurfaceProperty& sp1 = surface1->propertyAt(iy, iz);
	sp1.addBoundaryCondition(BoundaryCondition::kXWall);
	sp1.addBoundaryCondition(BoundaryCondition::kConstantTheta);
	sp1.setTemperature(10.0);
      }
      if (surface2) {
	SurfaceProperty& sp2 = surface2->propertyAt(iy, iz);
	sp2.addBoundaryCondition(BoundaryCondition::kXWall);
	sp2.addBoundaryCondition(BoundaryCondition::kConstantTheta);
	sp2.setTemperature(20.0);
      }
    }
  }
  surface1 = mEnvironment->surfaceState(BoundaryCell::kYLow);
  surface2 = mEnvironment->surfaceState(BoundaryCell::kYHigh);
  for (iz=0; iz<nz; iz++) {
    for (ix=0; ix<nx; ++ix) {
      if (surface1) {
	SurfaceProperty& sp1 = surface1->propertyAt(iz, ix);
	sp1.addBoundaryCondition(BoundaryCondition::kYWall);
      }
      if (surface2) {
	SurfaceProperty& sp2 = surface2->propertyAt(iz, ix);
	sp2.addBoundaryCondition(BoundaryCondition::kYWall);
      }
    }
  }
  surface1 = mEnvironment->surfaceState(BoundaryCell::kZLow);
  surface2 = mEnvironment->surfaceState(BoundaryCell::kZHigh);
  for (ix=0; ix<nx; ++ix) {
    for (iy=0; iy<ny; iy++) {
      if (surface1) {
	SurfaceProperty& sp1 = surface1->propertyAt(ix, iy);
	sp1.addBoundaryCondition(BoundaryCondition::kZWall);
      }
      if (surface2) {
	SurfaceProperty& sp2 = surface2->propertyAt(ix, iy);
	sp2.addBoundaryCondition(BoundaryCondition::kZWall);
      }
    }
  }
  
}

