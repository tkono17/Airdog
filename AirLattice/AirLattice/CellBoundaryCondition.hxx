#ifndef __CellBoundaryCondition_hxx__
#define __CellBoundaryCondition_hxx__
/*
  CellBoundaryCondition.hxx
*/

class CellBoundaryCondition {
public:
  enum Type {
    PeriodicX1, PeriodicX2, PeriodicY1, PeriodicY2, PeriodicZ1, PeriodicZ2, 
    WallX1, WallX2, WallY1, WallY2, WallZ1, WallZ2, 
  };

  CellBoundaryCondition();
  ~CellBoundaryCondition();

  void set(CellBoundaryCondition::Type bctype, bool value=true);

  bool isSet(CellBoundaryCondition::Type bctype) const;

protected:
  unsigned int mFlags;
};

#endif // __CellBoundaryCondition_hxx__
