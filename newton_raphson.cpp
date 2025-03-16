
#include <iostream>
#include "eigen-3.4.0/Eigen/Eigen"
#include "eigen-3.4.0/Eigen/Dense"
#include "eigen-3.4.0/Eigen/Sparse"
#include "eigen-3.4.0/Eigen/SparseLU"

using Eigen::VectorXd;
using Eigen::SparseMatrix;
template<typename Func, typename Jacobian>
void newton_raphson(Func F, Jacobian J, VectorXd& x0, double tol=1e-13, int max_nb_iterations=100){
    /* x_1 = x_0 - f(x_0)/f'(x_0) where we will use Jacobians instead */
    int iteration = 1;
    double residue = tol + 1;
    VectorXd x1;
    while (residue > tol && iteration < max_nb_iterations) {
        Eigen::MatrixXd J_dense = J(x0); 
        SparseMatrix<double> Jx0 = J_dense.sparseView();
        Eigen::SparseLU<SparseMatrix<double>> solver;
        solver.compute(Jx0);
        x1 = x0- solver.solve(F(x0));
        residue=(x1-x0).squaredNorm();
        x0 = x1;

        std::cout << "Residue after iteration " << iteration << " is " << residue << std::endl;
        iteration++;
    }
}