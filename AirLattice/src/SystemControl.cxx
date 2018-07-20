/*
  SystemControl.cxx
*/
#include <iostream>
#include "AirLattice/SystemControl.hxx"
#include "AirLattice/Constants.hxx"
#include "AirLattice/AirProperty.hxx"
#include "AirLattice/SurfaceState.hxx"

SystemControl::SystemControl() {
  mEnvironment = 0;
}

SystemControl::SystemControl(Environment* env, double dt, double ntimepoints) {
  mEnvironment = env;
  mDeltaTime = dt;
  mNTimePoints = ntimepoints;
}

SystemControl::~SystemControl() {
}

void SystemControl::updateEnvironment(float dt) {
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
	if (ix == 0) ix1 = ns[0] - 1;
	if (ix2 == ns[0]) ix2 = 0;
	if (iy == 0) iy1 = ns[1] - 1;
	if (iy2 == ns[1]) iy2 = 0;
	if (iz == 0) iz1 = 0;
	if (iz2 == ns[2]) iz2 = ns[2] - 1;

	mAP0 = &atmosphere.propertyAt(ix, iy, iz);
	mAPx1 = &atmosphere.propertyAt(ix1, iy, iz);
	mAPx2 = &atmosphere.propertyAt(ix2, iy, iz);
	mAPy1 = &atmosphere.propertyAt(ix, iy1, iz);
	mAPy2 = &atmosphere.propertyAt(ix, iy2, iz);
	mAPz1 = &atmosphere.propertyAt(ix, iy, iz1);
	mAPz2 = &atmosphere.propertyAt(ix, iy, iz2);

	updateTemperature();
	updatePressure();
	updateVelocity();
	updateDensity();
      }
    }
  }
}

void SystemControl::updateDensity() {
  double rho = mAP0->rho();

  double dudx=(mAPx2->u()-mAP0->u())/(mDeltaX);
  double dvdy=(mAPy2->v()-mAP0->v())/(mDeltaY);
  double dwdz=(mAPz2->w()-mAP0->w())/(mDeltaZ);

  rho -= rho*(dudx + dvdy + dwdz)*mDeltaTime;
  mAP0->setDensity(rho);
}

void SystemControl::updateVelocity() {
  int i;
  double v=0.0;
  double v0[3];
  double a=0.0;
  double dx=0.0;
  double dp=0.0;
  double fg=0.0;
  double rho = mAP0->rho();

  const double* vp = mAP0->velocity();
  for (i=0; i<2; ++i) {
    v0[i] = vp[i];
  }

  for (i=0; i<2; ++i) {
    v = v0[i];
    a = 0.0;

    a -= v0[0]*( mAPx2->velocity()[i] - vp[i] )/mDeltaX;
    a -= v0[1]*( mAPy2->velocity()[i] - vp[i] )/mDeltaY;
    a -= v0[2]*( mAPz2->velocity()[i] - vp[i] )/mDeltaZ;
    if (i==0) {
      dx = mDeltaX;
      dp = mAPx2->P() - mAP0->P();
      fg = 0.0;
    } else if (i==1) {
      dx = mDeltaY;
      dp = mAPz2->P() - mAP0->P();
      fg = 0.0;
    } else if (i==2) {
      dx = mDeltaZ;
      dp = mAPz2->P() - mAP0->P();
      fg = -Constants::g;
    }
    a -= dp/dx/rho;
    a += fg;

    v += a*mDeltaTime;    
  }
}

void SystemControl::updateTemperature() {
  double heat[3];
  double hc = mAP0->heatConductivity();
  //  double specheat = mAP0->specificHeat();
  double heatcapacity = mAP0->heatCapacity();
  //  double rho = mAP0->rho();
  double Q=0.0;
  double t = mAP0->T();

  heat[0] = hc*( mAPx2->T() - 2.0*t + mAPx1->T() );
  heat[1] = hc*( mAPy2->T() - 2.0*t + mAPy1->T() );
  if (mXYZIndex[2] != 0) {
    heat[2] = hc*( mAPz2->T() - 2.0*t + mAPz1->T() );
  } else {
    // This element touches the surface
    SurfaceState& surface = mEnvironment->surfaceState();
    SurfaceProperty& sp = surface.propertyAt(mXYZIndex[0], mXYZIndex[1]);
    heat[2] = hc*( mAPz2->T() - 2.0*t + sp.T());
  }
  Q = heat[0] + heat[1] + heat[2];

  t += Q*mDeltaTime/(heatcapacity);
  mAP0->setTemperature(t);
}

void SystemControl::updatePressure() {
  // PV = nRT <=> P = rho/M*RT
  double rho = mAP0->rho();
  double M = mAP0->molecularWeight();
  double t= mAP0->T();
  double p = rho/M*Constants::R*t;
  mAP0->setPressure(p);
}

int SystemControl::initialize() {
  const double* eSize = mEnvironment->elementSize();
  mAreaYZ = eSize[1]*eSize[2];
  mAreaZX = eSize[2]*eSize[0];
  mAreaXY = eSize[0]*eSize[1];
  mDeltaX = eSize[0];
  mDeltaY = eSize[1];
  mDeltaZ = eSize[2];
  return 0;
}

int SystemControl::run() {
  int it;
  double t;

  for (it=0; it<mNTimePoints; ++it) {
    std::cout << "System at time step " << it << std::endl;
    t = mDeltaTime*it;
    mTimeStep = it;

    std::cout << "  Updating environment" << std::endl;
    updateEnvironment(mDeltaTime);
    std::cout << "  Iteration done" << std::endl;
    if (it == 0 || it == (mNTimePoints-1) || 
	(it%mRecordInterval)==0) {
      mRecorder->save(it);
    }
  }

  return 0;
}
