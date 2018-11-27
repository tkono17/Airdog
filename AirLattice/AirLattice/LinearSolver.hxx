#ifndef __LinearSolver_hxx__
#define __LinearSolver_hxx__
/*
  LinearSolver.hxx
*/
#include "Eigen/Sparse"
#include "Eigen/Dense"

class LinearSolver {
public:
  LinearSolver();
  LinearSolver(int ndim);

  void assignMatrix(int nrow, int ncol, double value=0.0);
  void setMatrixElement(int i, int j, double value);

  void assignVector(int ndim, double value=0.0);
  void setVectorElement(int i, double value);

  void solve(const std::string& solver_name="BiCGSTAB");

  const Eigen::SparseMatrix<double>& matrixA() const { return mMatrixA; }
  const Eigen::VectorXd& vectorB() const { return mVectorB; }
  const Eigen::VectorXd& vectorX() const { return mVectorX; }
  Eigen::VectorXd& vectorX() { return mVectorX; }

  int ndim() const { return mN; }

protected:
  int mN;
  Eigen::SparseMatrix<double> mMatrixA;
  Eigen::VectorXd mVectorB;
  Eigen::VectorXd mVectorX;
};

#endif // __LinearSolver_hxx__
