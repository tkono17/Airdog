/*
  EnvRecorder.cxx
*/
#include <cstdio>
#include "AirLattice/EnvRecorder.hxx"
#include "TH2.h"

EnvRecorder::EnvRecorder() {
  mFilename = "envrec.root";
  mOutFile = 0;
  mSetName = "Set1";
  mEnvironment = 0;
}

EnvRecorder::EnvRecorder(const Environment* env) {
  mFilename = "envrec.root";
  mOutFile = 0;
  mSetName = "Set1";
  mEnvironment = env;
}

EnvRecorder::~EnvRecorder() {
}

void EnvRecorder::openFile(const std::string& fname) {
  mFilename = fname;
  mOutFile = TFile::Open(mFilename.c_str(), "RECREATE");
}

void EnvRecorder::closeFile() {
  if (mOutFile) {
    mOutFile->Write();
    mOutFile->Close();
    delete mOutFile;
    mOutFile = 0;
  }
}

void EnvRecorder::save(int timestep) {
  if (mSetName == "Set1") {
    saveSet1(timestep);
  }
}

void EnvRecorder::saveXZ(int iy, int timestep, PropertyType pt) {
  int ix, iz;
  int nx, ny, nz;
  const int* np = mEnvironment->NPoints();
  double value=0.0;
  char hname[200];

  nx = np[0];
  ny = np[1];
  nz = np[2];

  std::sprintf(hname, "h_%s_XZ_x%04d_t%04d", 
	       propertyName(pt).c_str(), iy, timestep);
  TH2F* h = new TH2F(hname, "", nx, 0, nx, nz, 0, nz);
  h->SetDirectory(mOutFile);

  for (ix=0; ix<nx; ++ix) {
    for (iz=0; iz<nz; ++iz) {
      value = mEnvironment->propertyValue(ix, iy, iz, pt);
      h->SetBinContent(ix+1, iz+1, value);
    }
  }
  h->Write();
  mOutFile->Write();
  delete h;
}

void EnvRecorder::saveXY(int iz, int timestep, PropertyType pt) {
  int ix, iy;
  int nx, ny, nz;
  const int* np = mEnvironment->NPoints();
  double value=0.0;
  char hname[200];

  nx = np[0];
  ny = np[1];
  nz = np[2];

  std::sprintf(hname, "h_%s_XY_z%04d_t%04d", 
	       propertyName(pt).c_str(), iz, timestep);
  TH2F* h = new TH2F(hname, "", nx, 0, nx, ny, 0, ny);
  h->SetDirectory(mOutFile);

  for (ix=0; ix<nx; ++ix) {
    for (iy=0; iy<ny; ++iy) {
      value = mEnvironment->propertyValue(ix, iy, iz, pt);
      h->SetBinContent(ix+1, iy+1, value);
    }
  }
  h->Write();
  mOutFile->Write();
  delete h;
}

void EnvRecorder::saveYZ(int ix, int timestep, PropertyType pt) {
  int iy, iz;
  int nx, ny, nz;
  const int* np = mEnvironment->NPoints();
  double value=0.0;
  char hname[200];

  nx = np[0];
  ny = np[1];
  nz = np[2];

  std::sprintf(hname, "h_%s_YZ_x%04d_t%04d", 
	       propertyName(pt).c_str(), ix, timestep);
  TH2F* h = new TH2F(hname, "", ny, 0, ny, nz, 0, nz);
  h->SetDirectory(mOutFile);

  for (iy=0; iy<ny; ++iy) {
    for (iz=0; iz<nz; ++iz) {
      value = mEnvironment->propertyValue(ix, iy, iz, pt);
      h->SetBinContent(iy+1, iz+1, value);
    }
  }
  h->Write();
  mOutFile->Write();
  delete h;
}

void EnvRecorder::saveSet1(int timestep) {
  saveXY(0, timestep, kRho);
  saveXY(0, timestep, kP);
  saveXY(0, timestep, kT);
  saveXY(0, timestep, kU);
  saveXY(0, timestep, kV);
  saveXY(0, timestep, kW);
  saveXZ(0, timestep, kRho);
  saveXZ(0, timestep, kP);
  saveXZ(0, timestep, kT);
}

