#ifndef __EnvRecorder_hxx__
#define __EnvRecorder_hxx__
/*
  EnvRecorder.hxx
*/
#include <string>
#include "AirLattice/Environment.hxx"
#include "AirLattice/PropertyType.hxx"
#include "TFile.h"

class EnvRecorder {
public:
  EnvRecorder();
  EnvRecorder(const Environment* env);
  ~EnvRecorder();

  void openFile(const std::string& fname);
  void closeFile();

  void setEnvironment(const Environment* x) { mEnvironment = x; }

  void save(int timestep);

  void saveXZ(int iy, int timestep, PropertyType pt);
  void saveXY(int iz, int timestep, PropertyType pt);
  void saveYZ(int ix, int timestep, PropertyType pt);

  void saveSet1(int timestep);

private:
  std::string mFilename;
  TFile* mOutFile;
  std::string mSetName;

  const Environment* mEnvironment;
};

#endif // __EnvRecorder_hxx__
