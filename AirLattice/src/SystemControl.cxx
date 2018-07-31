/*
  SystemControl.cxx
*/
#include <cstdio>
#include <iostream>
#include <functional>
#include "AirLattice/SystemControl.hxx"
#include "AirLattice/Constants.hxx"
#include "AirLattice/AirProperty.hxx"
#include "AirLattice/SurfaceState.hxx"

SystemControl::SystemControl() {
  mEnvironment = 0;
  mDoPrint=false;
}

SystemControl::SystemControl(Environment* env, double dt, double ntimepoints) {
  mEnvironment = env;
  mDeltaTime = dt;
  mNTimePoints = ntimepoints;
  mDoPrint=false;
}

SystemControl::~SystemControl() {
}

void SystemControl::updateEnvironment(float dt) {
  loop(std::mem_fun(&SystemControl::updateTemperature) );
  loop(std::mem_fun(&SystemControl::updatePressure) );
  loop(std::mem_fun(&SystemControl::updateVelocity) );
  loop(std::mem_fun(&SystemControl::updateDensity) );
  loop(std::mem_fun(&SystemControl::updateVariables) );
}

void SystemControl::loop(std::mem_fun_t<void, SystemControl> action) {
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

	mXYZIndex[0] = ix;
	mXYZIndex[1] = iy;
	mXYZIndex[2] = iz;

	mAP0 = &atmosphere.propertyAt(ix, iy, iz);
	mAPx1 = &atmosphere.propertyAt(ix1, iy, iz);
	mAPx2 = &atmosphere.propertyAt(ix2, iy, iz);
	mAPy1 = &atmosphere.propertyAt(ix, iy1, iz);
	mAPy2 = &atmosphere.propertyAt(ix, iy2, iz);
	mAPz1 = &atmosphere.propertyAt(ix, iy, iz1);
	mAPz2 = &atmosphere.propertyAt(ix, iy, iz2);

	mDoPrint = false;
	if (ix==500 && iy==0 && iz==0) {
	  mDoPrint = true;
	  mAP0->print();
	}
	action(this);
      }
    }
  }
}

void SystemControl::updateDensity() {
  mAP0->setDensityFromPT();
  // double rho = mAP0->rho();

  // double dudx=(mAPx2->u()-mAPx1->u())/(2*mDeltaX);
  // double dvdy=(mAPy2->v()-mAPy1->v())/(2*mDeltaY);
  // double dwdz=(mAPz2->w()-mAPz1->w())/(2*mDeltaZ);

  // if (mDoPrint) {
  //   std::cout << "Update rho: rho=" << rho
  // 	      << " dudx=" << dudx << ", dvdy=" << dvdy << ", dwdz=" << dwdz
  // 	      << " dfrac=" << (dudx + dvdy + dwdz)*mDeltaTime
  // 	      << std::endl;
  // }
  // rho -= rho*(dudx + dvdy + dwdz)*mDeltaTime;
  // mAP0->setDensity(rho);
}

void SystemControl::updateVelocity() {
  int i;
  double v=0.0;
  double v0[3], v1[3];
  double a=0.0;
  double dx=0.0;
  double dp=0.0;
  double fg=0.0;
  double rho = mAP0->rhoSI();
  double p0 = mAP0->P();

  const double* vp = mAP0->velocity();
  for (i=0; i<3; ++i) {
    v0[i] = vp[i];
    v1[i] = vp[i];
  }

  for (i=0; i<3; ++i) {
    v = v0[i];
    a = 0.0;

    a -= v0[0]*( mAPx2->velocity()[i] - v0[i] )/mDeltaX;
    a -= v0[1]*( mAPy2->velocity()[i] - v0[i] )/mDeltaY;
    a -= v0[2]*( mAPz2->velocity()[i] - v0[i] )/mDeltaZ;
    if (i==0) {
      // dx = mDeltaX;
      // dp = mAPx2->P() - p0;
      dx = 2*mDeltaX;
      dp = mAPx2->P() - mAPx1->P();
      fg = 0.0;
    } else if (i==1) {
      // dx = mDeltaY;
      // dp = mAPy2->P() - p0;
      dx = 2*mDeltaY;
      dp = mAPy2->P() - mAPy1->P();
      fg = 0.0;
    } else if (i==2) {
      // dx = mDeltaZ;
      // dp = mAPz2->P() - p0;
      dx = 2*mDeltaZ;
      dp = mAPz2->P() - mAPz1->P();
      fg = -Constants::g;
    }
    a = -dp/(dx*rho);
    a += fg;

    if (mDoPrint) {
      std::cout << "dv[" << i << "]: dp=" 
		<< (-dp/dx/rho) << " fg=" << fg 
		<< " dp=" << dp << " dx=" << dx << " rho=" << rho
		<< std::endl;
    }
    v += a*mDeltaTime;
    v1[i] = v;
  }
  if (mDoPrint) {
     std::cout << " v1 = (" << v1[0] << ", " << v1[1] << ", " << v1[2] << ")"
	       << std::endl;
  }
  //  mAP0->setVelocity(v1[0], v1[1], v1[2]);
  mAP0->updateVariable(kU, v1[0]);
  mAP0->updateVariable(kV, v1[1]);
  mAP0->updateVariable(kW, v1[2]);
}

void SystemControl::updateTemperature() {
  double heat[3];
  double hc = mAP0->heatConductivity();
  double heatcapacity = mAP0->heatCapacity(elementVolume());
  double Q=0.0;
  double t = mAP0->T();

  heat[0] = hc*( mAPx2->T() - 2.0*t + mAPx1->T() )/(mDeltaX*mDeltaX);
  heat[1] = hc*( mAPy2->T() - 2.0*t + mAPy1->T() )/(mDeltaY*mDeltaY);
  if (mXYZIndex[2] != 0) {
    heat[2] = hc*( mAPz2->T() - 2.0*t + mAPz1->T() )/(mDeltaZ*mDeltaZ);
  } else {
    // This element touches the surface
    SurfaceState& surface = mEnvironment->surfaceState();
    SurfaceProperty& sp = surface.propertyAt(mXYZIndex[0], mXYZIndex[1]);
    heat[2] = hc*( mAPz2->T() - 2.0*t + sp.T())/(mDeltaZ*mDeltaZ);
  }
  Q = (heat[0] + heat[1] + heat[2])*elementVolume();

  t += Q*mDeltaTime/(heatcapacity);

  if (mDoPrint) {
    char s[200];
    std::cout << "Q=" << Q << ", dT = " << t-mAP0->T() << std::endl;
    std::cout << " hc=" << hc << ", dt=" << mDeltaTime 
	      << " cp=" << heatcapacity 
	      << ", dx=" << mDeltaX 
	      << ", dy=" << mDeltaY 
	      << ", dz=" << mDeltaZ 
	      << std::endl;
    std::sprintf(s, " dV=%f %f %e", 0.00002, elementVolume(), elementVolume());
    std::cout << s << std::endl;
  }
  // std::cout << "T(" <<mXYZIndex[0]<<","<<mXYZIndex[1]<<","<<mXYZIndex[2]<<")="
  // 	    << t << std::endl;
  //  mAP0->setTemperature(t);
  mAP0->updateVariable(kT, t);
}

void SystemControl::updatePressure() {
  // PV = nRT <=> P = rho/M*RT
  double n = mAP0->molDensitySI();
  double T= mAP0->T();
  double p = n*Constants::R*T;
  if (mDoPrint) {
    std::cout << "mol density: " << n << ", T=" << T << ", R=" << Constants::R
	      << ", molecularWeight=" << mAP0->molecularWeight()
	      << ", density=" << mAP0->density()
	      << std::endl;
    std::cout << "p: " << mAP0->P() << " -> " << p << std::endl;
  }
  mAP0->setPressure(p);
  mAP0->updateVariable(kP, p);
}

void SystemControl::updateVariables() {
  mAP0->updateVariables();
}

void SystemControl::applyInitialConditions() {
  return uniformHotspot1();

  // std::set<std::string>::const_iterator p1, p2, p3, p4, p5;

  // bool apply_air_Uniform=false;
  // bool apply_air_T25C=false;
  // bool apply_air_P1Atom=false;
  // bool apply_surface_T25C=false;
  // bool apply_surface_Hotspot30C=false;
  // p1 = mInitialCondition.find("AirUniform10km");
  // p2 = mInitialCondition.find("AirT25C");
  // p3 = mInitialCondition.find("AirP1Atom");
  // p4 = mInitialCondition.find("SurfaceaT25C");
  // p5 = mInitialCondition.find("SurfaceHotspot30C");
  // if (p1 != mInitialCondition.end() ) {
  //   apply_air_Uniform = true;
  // }
  // if (p2 != mInitialCondition.end() ) {
  //   apply_air_T25C = true;
  // }
  // if (p3 != mInitialCondition.end() ) {
  //   apply_air_P1Atom = true;
  // }
  // if (p4 != mInitialCondition.end() ) {
  //   apply_surface_T25C = true;
  // }
  // if (p5 != mInitialCondition.end() ) {
  //   apply_surface_Hotspot30C = true;
  // }


  // int ix, iy, iz;
  // const int* np = mEnvironment->NPoints();
  // Atmosphere& atmosphere = mEnvironment->atmosphere();
  // SurfaceState& surface = mEnvironment->surfaceState();
  // const double t0 = 273.15;

  // int nx = np[0];
  // int ny = np[1];
  // int nz = np[2];
  // int hsx = static_cast<int>(np[0]/2);
  // int hsy = static_cast<int>(np[1]/2);

  // for (ix=0; ix<nx; ++ix) {
  //   for (iy=0; iy<ny; ++iy) {
  //     for (iz=0; ix<nz; ++iz) {
  // 	AirProperty& ap = atmosphere.propertyAt(ix, iy, iz);
  // 	if (apply_air_uniform) {
  // 	  if (apply_air_T25C) {
  // 	    ap.setTemperature(25 + t0);
  // 	  }
  // 	  if (apply_air_P1Atom) {
  // 	    ap.setPressure(1013.0E+2);
  // 	  }
  // 	  ap.setVelocity(0.0, 0.0, 0.0);
  // 	  ap.setVaporPressure(0.0);
  // 	  ap.setWaterDensity(0.0);
  // 	}
  //     }
  //     // Set surface property (ix, iy)
  //     SufaceProperty& sp = surface.propertyAt(ix, iy);
  //     if (apply_surface_T25C) {
  // 	sp.setTemperature(25 + t0);
  //     }
  //     if (apply_surface_Hotspot30C) {
  // 	if ( (ix == hsx) && (iy == hsy) ) {
  // 	  sp.setTemperature(30 + t0);
  // 	}
  //     }
  //   }
  // }

}

void SystemControl::applyBoundaryConditions() {
  // std::set<std::string>::const_iterator p1, p2, p3;

  // p1 = mBoundaryCondition.find("PeriodicXY");
  // p2 = mBoundaryCondition.find("SurfaceT");
  // p3 = mBoundaryCondition.find("SurfaceHotspot");

}

int SystemControl::initialize() {
  std::cout << "Initializing the system control" << std::endl;

  const double* eSize = mEnvironment->elementSize();
  mAreaYZ = eSize[1]*eSize[2];
  mAreaZX = eSize[2]*eSize[0];
  mAreaXY = eSize[0]*eSize[1];
  mDeltaX = eSize[0];
  mDeltaY = eSize[1];
  mDeltaZ = eSize[2];

  std::cout << "  Apply initial conditions" << std::endl;
  applyInitialConditions();

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
    if (it == 0 || 
	(it%mRecordInterval)==0 ) {
      //    if (it == 0) {
      mRecorder->save(it);
    }
    updateEnvironment(mDeltaTime);
    std::cout << "  Iteration done" << std::endl;
    if (it == (mNTimePoints-1) ) {
      mRecorder->save(it);
    }
  }

  return 0;
}

void SystemControl::uniformHotspot1() {
  int ix, iy, iz;
  const int* np = mEnvironment->NPoints();
  Atmosphere& atmosphere = mEnvironment->atmosphere();
  SurfaceState& surface = mEnvironment->surfaceState();
  const double t0 = Constants::T0;

  std::cout << "UniformHotspot1 condition" << std::endl;

  int nx = np[0];
  int ny = np[1];
  int nz = np[2];
  int hsx = static_cast<int>(np[0]/2);
  int hsy = static_cast<int>(np[1]/2);
  double z=0.0;
  double dt = 200.0/nz;

  std::cout << "  Create a hotspot on the surface ("
	    << hsx << ", " << hsy << ")" << std::endl;

  for (ix=0; ix<nx; ++ix) {
    for (iy=0; iy<ny; ++iy) {
      for (iz=0; iz<nz; ++iz) {
	z = mEnvironment->z(iz);
	AirProperty& ap = atmosphere.propertyAt(ix, iy, iz);
	ap.setPressure(1013.0E+2);
	ap.setTemperature(25 + t0);
	ap.setVelocity(0.0, 0.0, 0.0);
	ap.setVaporPressure(0.0);
	ap.setWaterDensity(0.0);
	if (iz > nz/2) {
	  ap.setTemperature(25 + t0 - dt*(iz-nz/2) );
	  ap.setPressure(1.0E+2);
	}
	ap.setDensityFromPT();
      }

      // Set surface property (ix, iy)
      SurfaceProperty& sp = surface.propertyAt(ix, iy);
      //      sp.setTemperature(25 + t0);
      sp.setTemperature(100 + t0);
      if ( (ix == hsx) && (iy == hsy) ) {
	sp.setTemperature(30 + t0);
      }
    }
  }
  return;
}
