#ifndef __SurfaceProperty_hxx__
#define __SurfaceProperty_hxx__
/*
  SurfaceProperty.hxx
*/

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

  double T() const { return mTemperature; }

private:
  double mTemperature;
  SurfaceProperty::Type_t mType;
};

#endif // __SurfaceProperty_hxx__
