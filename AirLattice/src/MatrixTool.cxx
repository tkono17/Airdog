/*
  MatrixTool.cxx
*/
#include <iostream>
#include "AirLattice/MatrixTool.hxx"

int solveGausJordan(SparseMatrix& m, std::vector<double>& b) {
  int ok=0;
  int i, j, k;
  int nrows, ncols;
  double a, f, a0, a1;
  std::vector<int> columns0;
  std::vector<int> columns1;
  std::vector<int>::const_iterator p0, p1;

  nrows = m.NRows();
  ncols = m.NColumns();

  for (i=0; i<nrows; ++i) {
    SparseVector& rv0 = m.rowVector(i);
    a = rv0.getValue(i);
    if (a == 0.0) {
      // this may happen for obstacles in the system
      continue;
    }
    for (k=0; k<ncols; ++k) {
      a0 = m.getValue(i, k);
      m.setValue(i, k, a0/a);
    }
    b[i] /= a;

    columns0 = rv0.columnsWithValues();
    for (j=i+1; j<nrows; ++j) {
      SparseVector& rv1 = m.rowVector(j);
      if (rv1.hasValue(i) ) {
	f = rv1.getValue(i);
	for (p0=columns0.begin(); p0!=columns0.end(); ++p0) {
	  if ( (*p0) == i) {
	    rv1.setValue(i, 0.0);
	  } else {
	    a0 = rv0.getValue(*p0);
	    a1 = rv1.getValue(*p0);
	    rv1.setValue(*p0, a1 - f*a0);
	  }
	}
	b[j] -= f*b[i];
      }
    }
    for (j=0; j<i; ++j) {
      SparseVector& rv1 = m.rowVector(j);
      if (rv1.hasValue(i) ) {
	f = rv1.getValue(i)/a;
	for (p0=columns0.begin(); p0!=columns0.end(); ++p0) {
	  if ( (*p0) == i) {
	    rv1.setValue(i, 0.0);
	  } else {
	    a0 = rv0.getValue(*p0);
	    a1 = rv1.getValue(*p0);
	    rv1.setValue(*p0, a1 - f*a0);
	  }
	}
	b[j] -= f*b[i];
      }
    }
  }
  return ok;
}

int solveGausJordanBand3(SparseMatrix& m, std::vector<double>& b) {
  int ok=0;
  int i, j, k;
  int nrows, ncols;
  double a, f, ai, aj;
  std::vector<int> columns0;
  std::vector<int> columns1;
  std::vector<int>::const_iterator p0, p1;

  nrows = m.NRows();
  ncols = m.NColumns();

  for (i=0; i<nrows; ++i) {
    a = m.getValue(i, i);
    if (a == 0.0) continue;

    // Normalize row i
    for (k=0; k<ncols; ++k) {
      ai = m.getValue(i, k);
      m.setValue(i, k, ai/a);
    }
    b[i] /= a;

    // Subtract rowi*fji from rowj
    for (j=0; j<nrows; ++j) {
      if (j==i) continue;
      f = m.getValue(j, i);
      m.setValue(j, i, 0.0);
      for (k=i+1; k<(i+3); ++k) {
	if (k > (ncols-1) ) break;
	ai = m.getValue(i, k);
	aj = m.getValue(j, k);
	m.setValue(j, k, aj - f*ai);
      }
      b[j] -= f*b[i];
    }
  }
  return ok;
}

void printMatrix(const SparseMatrix& m) {
  int i, j;
  std::cout << "SparseMatrix " << m.NRows() << "x" << m.NColumns()
	    << std::endl;
  for (i=0; i<m.NRows(); ++i) {
    std::cout << "[ ";
    for (j=0; j<m.NColumns(); ++j) {
      std::cout << m.getValue(i, j) << " ";
    }
    std::cout << "]" << std::endl;
  }
}

void printVector(const std::vector<double>& v) {
  int i, j;
  std::cout << "vector " << v.size() << std::endl;
  std::cout << "[ ";
  for (i=0; i<v.size(); ++i) {
    std::cout << v[i] << " ";
  }
  std::cout << "]" << std::endl;
}
