/*
  PropertyType.cxx
*/
#include <map>
#include "AirLattice/PropertyType.hxx"

std::string propertyName(PropertyType pt) {
  static std::map<PropertyType, std::string> pm;
  if (pm.size() == 0) {
    pm[kRho] = "Rho";
    pm[kP] = "P";
    pm[kT] = "T";
    pm[kU] = "U";
    pm[kV] = "V";
    pm[kW] = "W";
    pm[kVaporPressure] = "VaporPressure";
    pm[kWaterDensity] = "WaterDensity";
    pm[kSurf_T] = "Surf_T";
    pm[kSurf_Material] = "Surf_Material";
    pm[kUnknown] = "Unknown";
  }
  return pm[pt];
}
