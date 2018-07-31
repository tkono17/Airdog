#ifndef __Atmosphere_hxx__
#define __Atmosphere_hxx__
/*
  Atmosphere.hxx
*/
#include "AirLattice/AirProperty.hxx"

class Atmosphere {
public:
  Atmosphere();
  ~Atmosphere();

  int allocateData(int nx, int ny, int nz);
  int freeData();

  const int* NPoints() const { return mNPoints; }

  AirProperty& propertyAt(int ix, int iy, int iz);
  const AirProperty& propertyAt(int ix, int iy, int iz) const;

  void dump();

private:
  int mNPoints[3];

  AirProperty*** mFieldData;

};

#endif // __Atomosphere_hxx__
