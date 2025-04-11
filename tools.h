#pragma once
#define EIGEN_USE_MKL_ALL

#include <vector> 
#include <iostream> 
#include <string>
#include <Eigen/Sparse>
namespace E = Eigen;
namespace PNPNS {
	typedef Eigen::SparseMatrix<double, Eigen::RowMajor, long long int> CSRMatrix;
	typedef Eigen::Triplet<double, long long int> ETrip;

	typedef struct {
		double x;
		double y;
		double z;
	} double3;

	typedef struct {
		int x;
		int y;
		int z;
	} int3;

	typedef struct {
		long long int surid;
		long long int type;
		long long int xleft;
		long long int yleft;
		long long int zleft;
		long long int xright;
		long long int yright;
		long long int zright;
	} int64_8;

	typedef struct {
		long long int idx;
		double x;
		double y;
		double z;
	} int64_double3;

	typedef struct {
		double diel;
		double dx0;
		double dx1;
		double dy0;
		double dy1;
		double dz0;
		double dz1;
	} double7;

	struct Atoms {
		std::vector<double> x;
		std::vector<double> y;
		std::vector<double> z;
		std::vector<double> charge;
		std::vector<double> radius;
	};
}
void printA(PNPNS::CSRMatrix& A, int i);
bool checkCSRMatrix(PNPNS::CSRMatrix& A);
double donutdistance(double r, double o, double off, double a, double b, double c);
void donutnormal(double r, double o, double off, double a, double b, double c, double* norm);
bool cubicwithdonut(double r, double o, double off, double x, double y, double z, double dx);
void readPQR(PNPNS::Atoms& atoms, std::string& fn, double probe, double mx, double my, double mz);
int cubicwithpore(double x, double y, double z, double dx, double r, double l, double debye);
int pointnanopore(double x, double y, double z, double r, double l);


