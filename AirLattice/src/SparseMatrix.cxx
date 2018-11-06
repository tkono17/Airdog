/*
  SparseMatrix.cxx
*/
#include <iostream>
#include "AirLattice/SparseMatrix.hxx"
#include "AirLattice/SparseVector.hxx"

SparseMatrix::SparseMatrix() {
}

SparseMatrix::~SparseMatrix() {
}

void SparseMatrix::initialize(int nrow, int ncol) {
  mNRows = nrow;
  mNColumns = ncol;
  SparseVector sv;
  mRowVectors.assign(mNRows, sv);
}

void SparseMatrix::setValue(int row, int col, double value) {
  mRowVectors[row].setValue(col, value);
}

double SparseMatrix::getValue(int row, int col) const {
  return mRowVectors[row].getValue(col);
}

bool SparseMatrix::hasValue(int row, int col) const {
  return mRowVectors[row].hasValue(col);
}

void SparseMatrix::clear() {
  std::vector<SparseVector>::iterator p;
  for (p=mRowVectors.begin(); p!=mRowVectors.end(); ++p) {
    p->clear();
  }
}

