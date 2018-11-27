#ifndef __SurfaceProperty_hxx__
#define __SurfaceProperty_hxx__
/*
  SurfaceProperty.hxx
*/
#include "AirLattice/BoundaryCondition.hxx"

class SurfaceProperty {
public:
  enum Type_t {
    kWater, 
    kGround, 
  };
public:
  SurfaceProperty();
  ~SurfaceProperty();

  void setTemperature(double x) { mTemperature = x; }

  void addBoundaryCondition(BoundaryCondition::Type bctype);
  bool hasBoundaryCondition(BoundaryCondition::Type bctype) const;
  BoundaryCondition boundaryCondition() const { return mBC; }

  double T() const { return mTemperature; }

  void print() const;

private:
  double mTemperature;
  SurfaceProperty::Type_t mType;
  BoundaryCondition mBC;
};

#endif // __SurfaceProperty_hxx__
