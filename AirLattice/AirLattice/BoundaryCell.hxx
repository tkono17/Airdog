#ifndef __BoundaryCell_hxx__
#define __BoundaryCell_hxx__
/*
  BoundaryCell.hxx
*/
#include "AirLattice/BoundaryCondition.hxx"
#include "AirLattice/SurfaceState.hxx"
#include "AirLattice/AirProperty.hxx"

class Environment;

class BoundaryCell {
public:
  enum Location {
    kXLow, kXHigh, 
    kYLow, kYHigh, 
    kZLow, kZHigh
  };

public:
  BoundaryCell();
  BoundaryCell(Environment* env, int ix, int iy, int iz);
  ~BoundaryCell();

  double P() const;
  double theta() const;
  double u() const;
  double v() const;
  double w() const;

  bool isXLow() const { return ( (mLocationFlags>>kXLow) & 0x1) == 1; }
  bool isXHigh() const { return ( (mLocationFlags>>kXHigh) & 0x1) == 1; }
  bool isYLow() const { return ( (mLocationFlags>>kYLow) & 0x1) == 1; }
  bool isYHigh() const { return ( (mLocationFlags>>kYHigh) & 0x1) == 1; }
  bool isZLow() const { return ( (mLocationFlags>>kZLow) & 0x1) == 1; }
  bool isZHigh() const { return ( (mLocationFlags>>kZHigh) & 0x1) == 1; }

protected:
  void setLocationFlags();
  const AirProperty& innerProperty() const;
  const AirProperty& periodicProperty() const;

protected:
  Environment* mEnv;
  int mIx;
  int mIy;
  int mIz;
  SurfaceState* mSurface;
  SurfaceProperty* mSurfaceProperty;
  //
  unsigned int mLocationFlags;
};

#endif // __BoundaryCell_hxx__
