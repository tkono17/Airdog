#ifndef __SparseMatrix_hxx__
#define __SparseMatrix_hxx__
/*
  SparseMatrix.hxx
*/
#include <vector>
#include "SparseVector.hxx"

class SparseMatrix {
public:
  SparseMatrix();
  ~SparseMatrix();

  void initialize(int row, int col);

  int NRows() const { return mNRows; }
  int NColumns() const { return mNColumns; }

  void setValue(int row, int col, double value);
  double getValue(int row, int col) const;
  bool hasValue(int row, int col) const;

  std::vector<SparseVector>& rowVectors() { return mRowVectors; }

  SparseVector& rowVector(int row) { return mRowVectors[row]; }

  void clear();

protected:
  int mNRows;
  int mNColumns;
  std::vector<SparseVector> mRowVectors;
};

#endif // __SparseMatrix_hxx__
