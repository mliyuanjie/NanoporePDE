// NanoporePDE.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include "mesh.h"
#include "tools.h"
#include "operator.h"
#include "solver.h"
#include <Eigen/IterativeLinearSolvers>

int main()
{
    std::string fn = "E:/code/NanoporePDE/1AO6.pqr";
    std::string fnout = "E:/code/NanoporePDE/mesh.bin";
    std::string fnout2 = "E:/code/NanoporePDE/point.bin";
    std::string fnout3 = "E:/code/NanoporePDE/matrix.bin";
    std::string fnout4 = "E:/code/NanoporePDE/result.bin";
    std::string fnout5 = "E:/code/NanoporePDE/model.vtk";
    OctaTree mesh = OctaTree(10, 7, 5, 0, 0, 0, 1024);
    mesh.SAS(fn, 1.4, 100, 300, 20);
    mesh.generateIndex();
    mesh.check();
    mesh.save(fnout, fnout2);
    mesh.toVTK(fnout5);
    mesh.info();
    
    PDEOperator pnp(&mesh);
    pnp.setPhysics(2000, 20, 1, -1, 1.96e-9, 2.03e-9, 997, 0.0089);
    long long int n_elem = mesh.index.size();
    std::vector<double> x(n_elem * 7, 0);
    for (long long int i = 0; i < n_elem; ++i) {
        x[i + 2 * n_elem] = 1;
        x[i + n_elem] = 1;
        x[i + 6 * n_elem] = 1;
    }
    std::vector<double> y(mesh.index.size() * 7, 0);
    pnp.ForwardPNPNS(x.data(), y.data());
    std::ofstream outfile(fnout4, std::ios::binary);
    pnp.righthand(y.data(), 0, 0.01, -0.01);
    outfile.write(reinterpret_cast<const char*>(y.data()), y.size() * sizeof(double));
    outfile.close();
    PNPNS::CSRMatrix J(mesh.index.size() * 3, mesh.index.size() * 3);
    pnp.JacobianPNPNS(J, x.data());
    J.makeCompressed();
    std::ofstream f;
    f.open(fnout3, std::ios::binary);
    int n_csr = J.nonZeros();
    f.write(reinterpret_cast<const char*>(&n_csr), sizeof(long long int));
    for (int k = 0; k < J.outerSize(); ++k) {
        for (PNPNS::CSRMatrix::InnerIterator it(J, k); it; ++it) {
            long long int row = it.row();
            long long int col = it.col();
            double value = it.value();
            f.write(reinterpret_cast<const char*>(&row), sizeof(long long int));
            f.write(reinterpret_cast<const char*>(&col), sizeof(long long int));
            f.write(reinterpret_cast<const char*>(&value), sizeof(double));
        }
    }
    std::cout << "write down" << std::endl;
    f.close();
    //Eigen::IncompleteLUT<double, long long int> ilu;
    //ilu.compute(J);
    //if (ilu.info() != Eigen::Success) {
    //    std::cerr << "ILU decomposition failed!" << std::endl;
    //    return -1;
    //}
    PNPNSSolver solver = PNPNSSolver();
    Eigen::VectorXd dx = Eigen::VectorXd::Zero(mesh.index.size() * 3);
    //Eigen::VectorXd dx = Eigen::VectorXd::Random(mesh.index.size() * 3);
    std::cout << J.nonZeros() << std::endl;
    solver.solve(J, dx.data(), y.data(), mesh.index.size() * 3, J.nonZeros());
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
