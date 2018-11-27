#ifndef __BoundaryCondition_hxx__
#define __BoundaryCondition_hxx__
/*
  BoundaryCondition.hxx
*/

class BoundaryCondition {
public:
  enum Type {
    kXPeriodic, kYPeriodic, kZPeriodic, 
    kXWall, kYWall, kZWall, 
    kConstantTheta
  };

public:
  BoundaryCondition();
  ~BoundaryCondition();

  unsigned int flags() const { return mFlags; }

  void addType(BoundaryCondition::Type bctype);
  bool hasType(BoundaryCondition::Type bctype) const;


protected:
  unsigned int mFlags;

};

#endif // __BoundaryCondition_hxx__
