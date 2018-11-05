/*
  test_diffusion_gs.cxx
*/
#include <iostream>
#include <cstdio>
#include "AirLattice/MatrixTool.hxx"
#include "AirLattice/SparseMatrix.hxx"
#include "TH1.h"
#include "TFile.h"

int main(int argc, char* argv[]) {
  const int n=10000; // number of time steps
  const int np=100; // number of lattice points
  const double dt=0.01;
  const double dx=0.001;
  const double Tl=20.0;
  const double Tr=20.0;

  TFile* f = TFile::Open("diffusion_gs.root", "RECREATE");
  TH1* hists[n+1];

  double kappa=0.0241; // J.K-1.m-1.s-1
  double cp=1006; // J.K-1.kg-1
  double rho=1.166; // kg.m-3
  double a = kappa/(cp*rho); // m2.s-1
  // Initial temperature T0=0


  double xi = a*dt/(dx*dx);
  std::cout << "xi = " << xi << std::endl;

  int i, j, i1, i2;

  SparseMatrix A;
  std::vector<double> b;
  std::cout << "Setup the matrix and vector np=" << np << std::endl;
  A.initialize(np, np);
  b.assign(np, 0.0);

  std::cout << "Setup the matrix to solve Poisson's equation" << std::endl;
  // Stationary case (d/dt->0)
  for (i=0; i<np; ++i) {
    i1 = i - 1;
    i2 = i + 1;
    if (i == 0) {
      A.setValue(i, i, -3.0*xi);
      A.setValue(i, i2, 1.0*xi);
      b[i] = -2.0*Tl*xi;
    } else if (i == (np-1) ) {
      A.setValue(i, i1, 1.0*xi);
      A.setValue(i, i, -3.0*xi);
      b[i] = -2.0*Tr*xi;
    } else {
      A.setValue(i, i1, 1.0*xi);
      A.setValue(i, i, -2.0*xi);
      A.setValue(i, i2, 1.0*xi);
      b[i] = 0.0;
    }
  }
  // std::cout << "Matrix A:" << std::endl;
  // printMatrix(A);
  // std::cout << "Vector b:" << std::endl;
  // printVector(b);

  // std::cout << "Solve Poisson's equation with Gaus-Jordan method" << std::endl;
  solveGausSeidel(A, b);
  // std::cout << "Print matrix and vector after solving the equation" << std::endl;
  // std::cout << "Matrix A:" << std::endl;
  // printMatrix(A);
  // std::cout << "Vector b:" << std::endl;
  // printVector(b);

  std::cout << "done" << std::endl;
  hists[0] = new TH1F("h0", "Stationary problem", np, 0, np);
  std::cout << "Fill histogram" << std::endl;
  for (i=0; i<np; ++i) {
    hists[0]->SetBinContent(i+1, b[i]);
  }

  // Non-stationary case
  int it;
  std::vector<double> tvec, dtheta;
  for (i=0; i<np; ++i) {
    tvec.push_back(0.0);
    dtheta.push_back(0.0);
  }
  char hname[100];

  // Time steps
  for (it=0; it<n; ++it) {
    std::sprintf(hname, "h%04d", it);
    hists[it] = new TH1F(hname, "", np, 0, np);

    for (i=0; i<np; ++i) {
      i1 = i - 1;
      i2 = i + 1;
      if (i == 0) {
	dtheta[i] = xi*(tvec[i2] -3*tvec[i] + 2*Tl);
      } else if (i == (np-1) ) {
	dtheta[i] = xi*(tvec[i1] -3*tvec[i] + 2*Tr);
      } else {
	dtheta[i] = xi*(tvec[i2] -2*tvec[i] + tvec[i1]);
      }
    }
    for (i=0; i<np; ++i) {
      tvec[i] += dtheta[i];
      hists[it]->SetBinContent(i+1, tvec[i]);
    }
  }

  f->Write();
  f->Close();

  return 0;
}
