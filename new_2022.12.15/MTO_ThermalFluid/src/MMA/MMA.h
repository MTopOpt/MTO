#ifndef MMA_H
#define MMA_H

#include "mpi.h"
#include <vector>
#include <cmath>
// #include "fvCFD.H"

class MMA {

  public:
	MMA(int NvarLocal, int m);

  void MMAsolver(std::vector<double> &xval,
                 std::vector<double> &dfdx,
                 std::vector<double> &g,
                 std::vector<std::vector<double>> &dgdx);
                 
public:
  int n;// number of design variables
  int m;// number of constraints
  int iter;//iteration

  // some parameters set by user
  double raa0=1.0e-5;
  double asyminit=0.5;
  double asymdec=0.7;
  double asyminc=1.2;
  double albefa=0.1;
  double epsimin;
  double xmamieps=1.0e-5;
  double movelimit=0.5;
  double RobustAsymptotesType=0;
  std::vector<double> a,c;

  private:
	std::vector<double> y;
	double z;
	std::vector<double> lam, mu, s;
	std::vector<double>  alpha, beta, p0, q0, pij, qij, b, grad, hess;

  std::vector<double> xmax, xmin, low, upp, xold1, xold2;



private:
  void Update(std::vector<double> &xval,
              std::vector<double> &dfdx,
              std::vector<double> &g,
              std::vector<std::vector<double>> &dgdx);

  void GenSub(std::vector<double> &xval, std::vector<double> &dfdx,
               std::vector<double> &g, std::vector<std::vector<double>> &dgdx);

	void SolveDIP(std::vector<double> &x);

	void XYZofLAMBDA(std::vector<double> &x);

	void DualGrad(std::vector<double> &x);
	void DualHess(std::vector<double> &x);
	void DualLineSearch();
	double DualResidual(std::vector<double> &x, double epsi);

	static void Factorize(double *K, int n);
	static void Solve(double *K, double *x, int n);
};
#endif