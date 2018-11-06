#ifndef __SparseVector_hxx__
#define __SparseVector_hxx__
/*
  SparseVector.hxx
*/
#include <map>
#include <vector>

class SparseVector {
public:
  SparseVector();
  ~SparseVector();

  void setValue(int col, double value);
  double getValue(int col) const;
  bool isNonZero(int col) const;
  bool hasValue(int col) const;

  std::vector<int> columnsWithValues() const;

  void clear();

protected:
  std::map<int, double> mElements;
};

#endif // __SparseVector_hxx__
