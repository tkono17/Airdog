#ifndef __PropertyType_hxx__
#define __PropertyType_hxx__
/*
  PropertyType.hxx
*/
#include <string>

enum PropertyType {
  kRho, 
  kP, 
  kT, 
  kU, kV, kW, 
  kVaporPressure, 
  kWaterDensity, 
  kSurf_T, 
  kSurf_Material, 
  kUnknown
};

std::string propertyName(PropertyType pt);


#endif // __PropertyType_hxx__

