#pragma once
#include <iostream>
#include "mesh.h"
#include "tools.h"
#include <list>

struct PartialNode {
	long long int index;
	double value;
	PartialNode(long int i, double v) : index(i), value(v) {}
};

typedef std::list<PartialNode> PartialList;

class PDEOperator {
public:
	PDEOperator(OctaTree* mesh);
	OctaTree* mesh;
	int graph[6][4] = {
		{ 1, 3, 5, 7 },
		{ 0, 2, 4, 6 },
		{ 2, 3, 6, 7 },
		{ 0, 1, 4, 5 },
		{ 4, 5, 6, 7 },
		{ 0, 1, 2, 3 }
		
	};
	PartialList* JMatrix = nullptr;
	void JacobianPNPNS(PNPNS::CSRMatrix& J, double* x);
	void gradJ(TreeNode* node, PartialList& y, double* coeff);
	void gradJ2(TreeNode* node, PartialList& y, double* coeff);
	void laplJ(TreeNode* node, PartialList& y, double* coeff);
	void gradJ_direction(TreeNode* node, PartialList& y, double* direction);

	void ForwardPNPNS(double* x, double* y);
	void grad(TreeNode* node, double* x, double* y);
	void lapl(TreeNode* node, double* x, double* y, double* coeff);
	void grad_direction(TreeNode* node, double* x, double* y, double* direction);
	void setPhysics(double conc, double temperature, int z0, int z1, double Diff0, double Diff1, double rho, double viscosity);
	void righthand(double* x, double v0, double v1, double chargedensity);

private:
	double ee = 1.602176634e-19;
	double eps0 = 8.854187817e-12 * 78.5;
	double diep = 3.0 / 78.5;
	double kb = 1.380649e-23;
	double avogadro = 6.02214076e23;
	double inv_eps0 = 1.0 / eps0;
	double kbt = 1.380649e-23 * (273.15 + 25);
	double ekbt = ee / kbt;
	double c0 = 2000;
	double z0 = 1;
	double z1 = -1;
	double diffu0 = 1.96e-9;
	double diffu1 = 2.03e-9;
	double rho = 997;
	double viscosity = 0.0089;
	double P0 = 101325;

	//poisson source scale 
	double charge_coeff0 = 0;
	double charge_coeff1 = 0;
	//nernst planck couple poisson
	//this is z0 z1, 
	//nernst planck couple fluid
	double np_ns_coeff0 = 0;
	double np_ns_coeff1 = 0;
	//navier-stock couple with pnp
	double ns_pnp_coeff = 0;
	double Re_inv = 0;
};
