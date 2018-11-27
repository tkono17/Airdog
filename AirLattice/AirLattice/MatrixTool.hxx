#ifndef __MatrixTool_hxx__
#define __MatrixTool_hxx__
/*
  MatrixTool.hxx
*/
#include <vector>
#include "AirLattice/SparseMatrix.hxx"

/**
   Solve the matrix equation m*x = b by Gaus-Jordan method.
 */
int solveGausJordan(SparseMatrix& m, std::vector<double>& b);

/**
   Solve the matrix equation m*x = b by Gaus-Jordan method for 
   band matrix with width=3.
 */
int solveGausJordanBand3(SparseMatrix& m, std::vector<double>& b);

/**
   Solve the matrix equation m*x = b by SOR method
 */
int solveSOR(SparseMatrix& m, std::vector<double>& b, 
	     const std::vector<double>& xinit, double omega=1.2);

/**
   Solve the matrix equation m*x = b by Gaus-Seidel method.
 */
int solveGausSeidel(SparseMatrix& m, std::vector<double>& b);

void printMatrix(const SparseMatrix& m);

void printVector(const std::vector<double>& v);

#endif // __MatrixTool_hxx__
