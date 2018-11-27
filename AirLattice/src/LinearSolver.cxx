/*
  LinearSolver.cxx
*/
#include "AirLattice/LinearSolver.hxx"
#include "Eigen/IterativeLinearSolvers"

LinearSolver::LinearSolver() {
}

LinearSolver::LinearSolver(int ndim) {
  assignMatrix(ndim, ndim);
  assignVector(ndim);
}

void LinearSolver::assignMatrix(int nrow, int ncol, double value) {
  mN = nrow;
  mMatrixA = Eigen::SparseMatrix<double>(nrow, ncol);
  mMatrixA.setZero();
}

void LinearSolver::setMatrixElement(int i, int j, double value) {
  mMatrixA.insert(i, j) = value;
}

void LinearSolver::assignVector(int ndim, double value) {
  mVectorB = Eigen::VectorXd(ndim);
  mVectorB.setZero();
}

void LinearSolver::setVectorElement(int i, double value) {
  mVectorB(i) = value;
}

void LinearSolver::solve(const std::string& solver_name) {
  if (solver_name == "BiCGSTAB") {
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > BiCGSTAB;
    BiCGSTAB.compute(mMatrixA);
    mVectorX = BiCGSTAB.solve(mVectorB);
  } else if (solver_name == "ConjugateGradient") {
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > solver;
    solver.compute(mMatrixA);
    mVectorX = solver.solve(mVectorB);
  } else {
    mVectorX.setZero();
  }
}



