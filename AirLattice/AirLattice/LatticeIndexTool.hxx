#ifndef __LatticeIndexTool_hxx__
#define __LatticeIndexTool_hxx__
/*
  LatticeIndexTool.hxx
*/
#include <vector>
#include "AirLattice/Environment.hxx"

class LatticeIndexTool {
public:
  LatticeIndexTool();
  ~LatticeIndexTool();

  int initialize(Environment* env);

  int Nx() const { return mNx; }
  int Ny() const { return mNy; }
  int Nz() const { return mNz; }
  int NCells() const { return mNCells; }
  
  int indexOfCell(int ix, int iy, int iz) const;

  std::vector<int> cellIndexAt(int icell) const;

protected:
  Environment* mEnvironment;

  // Auxiliary variables
  int mNx;
  int mNy;
  int mNz;
  int mNCells; // Total number of parameters (naively Nx*Ny*Nz)
};

#endif // __LatticeIndexTool_hxx__
