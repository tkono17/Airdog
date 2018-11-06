/*
  SparseVector.cxx
*/
#include "AirLattice/SparseVector.hxx"

SparseVector::SparseVector() {
}

SparseVector::~SparseVector() {
}

void SparseVector::setValue(int col, double value) {
  mElements[col] = value;
}

double SparseVector::getValue(int col) const {
  std::map<int, double>::const_iterator p=mElements.find(col);
  if (p == mElements.end() ) {
    return 0.0;
  } else {
    return (p->second);
  }
}

bool SparseVector::isNonZero(int col) const {
  bool x=false;
  std::map<int, double>::const_iterator p;
  if ( (p=mElements.find(col) ) == mElements.end() ) {
    x = false;
  } else {
    if ( (p->second) != 0.0) {
      x = true;
    }
  }
  return x;
}

bool SparseVector::hasValue(int col) const {
  bool x=false;
  std::map<int, double>::const_iterator p;
  if ( (p=mElements.find(col) ) == mElements.end() ) {
    x = false;
  } else {
    x = true;
  }
  return x;
}

std::vector<int> SparseVector::columnsWithValues() const {
  std::vector<int> v;
  std::map<int, double>::const_iterator p;
  for (p=mElements.begin(); p!=mElements.end(); ++p) {
    v.push_back(p->first);
  }
  return v;
}

void SparseVector::clear() {
  mElements.clear();
}
