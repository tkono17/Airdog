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
  loop(std::mem_fun(&SystemControl::updateVelocity) );
  loop(std::mem_fun(&SystemControl::updateDensity) );
  loop(std::mem_fun(&SystemControl::updatePressure) );
  loop(std::mem_fun(&SystemControl::updateVariables) );

  Atmosphere& atmosphere = mEnvironment->atmosphere();
  std::cout << "Field at (ix, iy, iz)=(49-51, 0, 1)" << std::endl;
  atmosphere.propertyAt(49, 0, 1).print();
  atmosphere.propertyAt(50, 0, 1).print();
  atmosphere.propertyAt(51, 0, 1).print();
  std::cout << "Field at (ix, iy, iz)=(49-51, 0, 0)" << std::endl;
  atmosphere.propertyAt(49, 0, 0).print();
  atmosphere.propertyAt(50, 0, 0).print();
  atmosphere.propertyAt(51, 0, 0).print();
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
	if (iz == 0) iz1 = -1;
	if (iz2 == ns[2]) iz2 = -1;

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
	mAPx1 = &atmosphere.propertyAt(ix1, iy, iz);
	mAPx2 = &atmosphere.propertyAt(ix2, iy, iz);
	mAPy1 = &atmosphere.propertyAt(ix, iy1, iz);
	mAPy2 = &atmosphere.propertyAt(ix, iy2, iz);
	if (iz1 < 0) {
	  mAPz1 = 0;
	} else {
	  mAPz1 = &atmosphere.propertyAt(ix, iy, iz1);
	}
	if (iz2 < 0) {
	  mAPz2 = 0;
	} else {
	  mAPz2 = &atmosphere.propertyAt(ix, iy, iz2);
	}

	action(this);
      }
    }
  }
}

void SystemControl::updateDensity() {
  //mAP0->setDensityFromPT();
  double rho = mAP0->rho();

  double drhou=0.0;
  double drhov=0.0;
  double drhow=0.0;
  drhou = -d(std::mem_fun<double, AirProperty>(&AirProperty::rhoU), 
	     mAP0, mAPx1, mAPx2, mDeltaX, 0.0, 0.0);
  drhov = -d(std::mem_fun<double, AirProperty>(&AirProperty::rhoV),
	     mAP0, mAPy1, mAPy2, mDeltaY, 0.0, 0.0);
  drhow = -d(std::mem_fun<double, AirProperty>(&AirProperty::rhoW),
	     mAP0, mAPz1, mAPz2, mDeltaZ, 0.0, 0.0);
  if (mDoPrint) {
    std::cout << "Update rho: rho=" << rho
   	      << " drhou=" << drhou
   	      << ", drhov=" << drhov
   	      << ", drhow=" << drhow
   	      << " drho=" << (drhou+drhov+drhow)*mDeltaTime
   	      << std::endl;
  }
  rho += (drhou+drhov+drhow)*mDeltaTime;
  //  mAP0->setDensity(rho);
  mAP0->updateVariable(kRho, rho);
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
  double advec=0.0;

  const double* vp = mAP0->velocity();
  for (i=0; i<3; ++i) {
    v0[i] = vp[i];
    v1[i] = vp[i];
  }

  for (i=0; i<3; ++i) {
    v = v0[i];
    a = 0.0;

    // if (i == 0) {
    //   advec = -vd(std::mem_fun<double, AirProperty>(&AirProperty::u), 0.0, 0.0);
    // } else if (i == 1) {
    //   advec = -vd(std::mem_fun<double, AirProperty>(&AirProperty::v), 0.0, 0.0);
    // } else if (i == 2) {
    //   advec = -vd(std::mem_fun<double, AirProperty>(&AirProperty::w), 0.0, 0.0);
    // }
    a -= v0[0]*( mAPx2->velocity()[i] - mAPx1->velocity()[i])/mDeltaX/2;
    a -= v0[1]*( mAPy2->velocity()[i] - mAPy1->velocity()[i])/mDeltaY/2;
    if (isAtBoundarySite(kZl) ) {
      a = 0.0;
    } else if (isAtBoundarySite(kZh) ) {
      a -= v0[2]*( mAP0->velocity()[i] - mAPz1->velocity()[i])/mDeltaZ/2;
    } else {
      a -= v0[2]*( mAPz2->velocity()[i] - mAPz1->velocity()[i])/mDeltaZ/2;
    }

    if (i==0) {
      // dx = mDeltaX;
      // dp = mAPx2->P() - p0;
      dp = -d_d(std::mem_fun<double, AirProperty>(&AirProperty::P), 
		mAP0, mAPx1, mAPx2, mDeltaX)/rho;
      fg = 0.0;
    } else if (i==1) {
      // dx = mDeltaY;
      // dp = mAPy2->P() - p0;
      dp = -d_d(std::mem_fun<double, AirProperty>(&AirProperty::P), 
		mAP0, mAPy1, mAPy2, mDeltaY)/rho;
      fg = 0.0;
    } else if (i==2) {
      // dx = mDeltaZ;
      // dp = mAPz2->P() - p0;
      if (isAtInnerSite()) {
	dp = -d_d(std::mem_fun<double, AirProperty>(&AirProperty::P), 
		  mAP0, mAPz1, mAPz2, mDeltaZ, 0.0, 0.0)/rho;
      } else if (isAtBoundarySite(kZl)) {
	// dp = -dForward(std::mem_fun<double, AirProperty>(&AirProperty::P), 
	// 	       mAP0, mAPz1, mAPz2, mDeltaZ, 0.0, 0.0)/rho;
	dp = -(mAP0->w() - 0.0)/mDeltaZ/rho;
      } else if (isAtBoundarySite(kZh)) {
	dp = -d_d(std::mem_fun<double, AirProperty>(&AirProperty::P), 
		   mAP0, mAPz1, mAPz2, mDeltaZ, 0.0, 0.0)/rho;
      }
      fg = -Constants::g;
    }
    a += advec;
    a += dp;
    a += fg;
    if (mDoPrint) {
      std::cout << "dv[" << i << "]: rho(SI)=" << rho << " dp=" 
		<< dp << " advection=" << advec << " fg=" << fg 
		<< " dx=" << dx << " rho=" << rho
		<< std::endl;
    }
    if (isAtBoundarySite(kZl) && i==2) {
      v = 0.0;
    } else {
      v += a*mDeltaTime;
    }
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

  SurfaceState& surface = mEnvironment->surfaceState();
  SurfaceProperty& sp = surface.propertyAt(mXYZIndex[0], mXYZIndex[1]);
  double bzl = sp.T();
  double bzh = mAP0->T();
  double advec=0.0;

  Q = laplacian(std::mem_fun<double, AirProperty>(&AirProperty::T), bzl, bzh);
  advec = -vd(std::mem_fun<double, AirProperty>(&AirProperty::T), bzl, bzh);

  Q = hc*Q*elementVolume();
  t += (Q/heatcapacity + advec)*mDeltaTime;

  if (mDoPrint) {
    char s[200];
    std::cout << " HeatConductivity(SI)=" << hc << ", dt=" << mDeltaTime 
	      << " HeatCapacity(SI)=" << heatcapacity 
	      << " advection=" << advec
	      << " Heat(in)=" << Q
	      << " dT=" << (t - mAP0->T() )
	      << ", dx=" << mDeltaX 
	      << ", dy=" << mDeltaY 
	      << ", dz=" << mDeltaZ 
	      << ", mAPz1=" << mAPz1 << ", surfaceT=" << bzl
	      << std::endl;
    std::sprintf(s, " VolumeElement=%f %f %e", 0.00002, elementVolume(), elementVolume());
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

  mDoPrint = false;
  for (it=0; it<mNTimePoints; ++it) {
    std::cout << "System at time step " << it << std::endl;
    t = mDeltaTime*it;
    mTimeStep = it;

    //    std::cout << "  Updating environment" << std::endl;
    if (it == 0 || 
	(it%mRecordInterval)==0 ) {
      //    if (it == 0) {
      mRecorder->save(it);
    }

    updateEnvironment(mDeltaTime);
    //    std::cout << "  Iteration done" << std::endl;
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
	// if (iz > nz/2) {
	//   ap.setTemperature(25 + t0 - dt*(iz-nz/2) );
	//   ap.setPressure(1.0E+2);
	// }
	ap.setDensityFromPT();
      }

      // Set surface property (ix, iy)
      SurfaceProperty& sp = surface.propertyAt(ix, iy);
      //      sp.setTemperature(25 + t0);
      sp.setTemperature(25 + t0);
      // if ( (ix == hsx) && (iy == hsy) ) {
      // 	sp.setTemperature(100 + t0);
      // }
    }
  }
  return;
}

double SystemControl::d(std::const_mem_fun_t<double, AirProperty> var, 
			AirProperty* s0, 
			AirProperty* sm, AirProperty* sp, 
			double dx, 
			double mBoundaryValue, double pBoundaryValue) {
  double df=0.0;

  //std::cout << "s0 " << s0 << " sm " << sm << " sp " << sp << std::endl;
  if (s0 && sm && sp) {
    //df = (var(sp) - var(sm))/(2.0*dx);
    df = (var(sp) - var(s0))/dx;
  } else if (s0 && (sm==0) && sp) {
    df = ( (var(sp)+var(s0) )/2.0 - mBoundaryValue)/dx;
  } else if (s0 && sm && (sp==0) ) {
    df = (pBoundaryValue - (var(s0)+var(sm) )/2.0)/dx;
  }
  return df;
}

double SystemControl::dForward(std::const_mem_fun_t<double, AirProperty> var, 
			       AirProperty* s0, 
			       AirProperty* sm, AirProperty* sp, 
			       double dx, 
			       double mBoundaryValue, double pBoundaryValue) {
  double df=0.0;

  if (s0 && sp) {
    df = (var(sp) - var(s0))/dx;
  }
  return df;
}

double SystemControl::dBackward(std::const_mem_fun_t<double, AirProperty> var, 
				AirProperty* s0, 
				AirProperty* sm, AirProperty* sp, 
				double dx, 
				double mBoundaryValue, double pBoundaryValue) {
  double df=0.0;

  if (s0 && sm) {
    df = (var(s0) - var(sm))/dx;
  }
  return df;
}

double SystemControl::d2(std::const_mem_fun_t<double, AirProperty> var, 
			 AirProperty* s0, 
			 AirProperty* sm, AirProperty* sp, 
			 double dx, 
			 double mBoundaryValue, double pBoundaryValue) {
  double df=0.0;

  if (s0 && sm && sp) {
    df = (var(sp) -2*var(s0) + var(sm))/(dx*dx);
  } else if (s0 && (sm==0) && sp) {
    df = (var(sp) -3*var(s0) + 2*mBoundaryValue)/(dx*dx);
  } else if (s0 && sm && (sp==0) ) {
    df = (2*pBoundaryValue - 3*var(s0) + var(sm) )/(dx*dx);
  }
  return df;
}

double SystemControl::vd(std::const_mem_fun_t<double, AirProperty> var, 
			 AirProperty* s0, 
			 AirProperty* sm, AirProperty* sp, 
			 double v, double dx, 
			 double mBoundaryValue, double pBoundaryValue) {
  double df=0.0;
  if ( (sm==0) || (sp==0) ) {
    df = v*d(var, s0, sm, sp, dx, mBoundaryValue, pBoundaryValue);
  } else {
    if (v > 0.0) {
      df = v*d(var, s0, sm, 0, dx, mBoundaryValue, pBoundaryValue);
    } else {
      df = v*d(var, s0, 0, sp, dx, mBoundaryValue, pBoundaryValue);
    }
  }
  return df;
}

double SystemControl::vd(std::const_mem_fun_t<double, AirProperty> var, 
			 double bzl, double bzh) {
  double y=0.0;

  y += vd(var, mAP0, mAPx1, mAPx2, mAP0->u(), mDeltaX);
  y += vd(var, mAP0, mAPy1, mAPy2, mAP0->v(), mDeltaY);
  y += vd(var, mAP0, mAPz1, mAPz2, mAP0->w(), mDeltaZ, bzl, bzh);
  return y;
}

double SystemControl::d_d(std::const_mem_fun_t<double, AirProperty> var, 
			  AirProperty* s0, AirProperty* sm, AirProperty* sp, 
			  double dx, 
			  double mBoundaryValue, double pBoundaryValue) {
  double y=0.0;
  y = (var(s0)-var(sm) )/dx;
  return y;
}

double SystemControl::vd_d(std::const_mem_fun_t<double, AirProperty> var, 
			   AirProperty* s0, AirProperty* sm, AirProperty* sp, 
			   double v, double dx, 
			   double mBoundaryValue, double pBoundaryValue) {
  double df=0.0;
  double f1=0.0;
  double f2=0.0;
  if (sm == 0) {
    df = v*(pBoundaryValue - var(s0) )/dx;
  } else if (sp == 0) {
    df = v*(var(s0) - mBoundaryValue)/dx;
  } else {
    f1 = (var(sm) + var(s0) )/2.0;
    f2 = (var(s0) + var(sp) )/2.0;
    df = v*(f2 - f1)/dx;
  }
  return df;
}

double SystemControl::vd_d(std::const_mem_fun_t<double, AirProperty> var, 
			   double bzl, double bzh) {
  double y=0.0;
  y += vd_d(var, mAP0, mAPx1, mAPx2, mAP0->u(), mDeltaX);
  y += vd_d(var, mAP0, mAPy1, mAPy2, mAP0->v(), mDeltaY);
  y += vd_d(var, mAP0, mAPz1, mAPz2, mAP0->w(), mDeltaZ, bzl, bzh);
  return y;
}

double SystemControl::laplacian(std::const_mem_fun_t<double, AirProperty> var, 
				double bzl, double bzh) {
  double y=0.0;
  y += d2(var, mAP0, mAPx1, mAPx2, mDeltaX);
  y += d2(var, mAP0, mAPy1, mAPy2, mDeltaY);
  y += d2(var, mAP0, mAPz1, mAPz2, mDeltaZ, bzl, bzh);
  return y;
}

bool SystemControl::isAtInnerSite() const {
  return (mAPx1 && mAPx2 && mAPy1 && mAPy2 && mAPz1 && mAPz2);
}

bool SystemControl::isAtBoundarySite(SystemControl::boundary_t bt) const {
  switch (bt) {
  case kXl:
    return (mAPx1 == 0);
    break;
  case kXh:
    return (mAPx2 == 0);
    break;
  case kYl:
    return (mAPy1 == 0);
    break;
  case kYh:
    return (mAPy2 == 0);
    break;
  case kZl:
    return (mAPz1 == 0);
    break;
  case kZh:
    return (mAPz2 == 0);
    break;
  default:
    return false;
  }
  return false;
}

