/*
  test_cellindex.cxx
*/
#include <iostream>
#include "AirLattice/LatticeIndexTool.hxx"

int main(int argc, char* argv) {
  LatticeIndexTool tool;
  Environment env;
  int nx=100;
  int ny=2;
  int nz=30;

  env.setSystemSize(10.0, 1.0, 10.0);
  env.setNPoints(nx, ny, nz);
  env.initialize();

  tool.initialize(&env);

  int ix, iy, iz;

  std::cout << "(Nx, Ny, Nz) = (" << nx << ", " << ny << ", " << nz << ")"
	    << std::endl;
  std::cout << "Total number of cells: " << tool.NCells() << std::endl;

  int naccum=0;
  int n_by_tool=0;
  bool ok=true;

  for (ix=0; ix<nx; ++ix) {
    for (iy=0; iy<ny; ++iy) {
      for (iz=0; iz<nz; ++iz) {
	n_by_tool = tool.indexOfCell(ix, iy, iz);
	if (n_by_tool != naccum) {
	  std::cout << "Cell number is not as expected, ("
		    << ix << ", " << iy << ", " << iz << ") -> "
		    << n_by_tool << "(expected=" << naccum << std::endl;
	  ok = false;
	}
	std::vector<int> v = tool.cellIndexAt(n_by_tool);
	if (n_by_tool != naccum) {
	  std::cout << "Cell number is not as expected, "
		    << n_by_tool << " -> ("
		    << v[0] << ", " << v[1] << ", " << v[2] 
		    << ") [expected=("
		    << ix << ", " << iy << ", " << iz << ")]"
		    << std::endl;
	  ok = false;
	}
	naccum ++;
      }
    }
  }
  if (ok) {
    std::cout << "LatticeIndexTool is working as expected" << std::endl;
  }
  return 0;
}
