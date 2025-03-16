#include <iostream>
#include <cmath>
#include "../matrix.cpp" 
#include "../read_file.hpp"
using namespace Eigen;
using namespace std;
/*----------------------------------------------------------------------------------------
This is a file for testing the correctness of the assembling of our matrix for the system
-----------------------------------------------------------------------------------------*/
void matrix_comparison(const MatrixXd& matrix_computed, const MatrixXd& matrix_expected, double tolerance) {

    if ((matrix_computed.rows()!=matrix_expected.rows()) || (matrix_computed.cols()!=matrix_expected.cols())) {
        throw runtime_error("Sorry the rows and columns of the matrices are not equal, there is an error somewhere");
    }
    
    for (int i = 0; i<matrix_computed.rows(); ++i) {
        for (int j = 0; j< matrix_computed.cols(); ++j) {
            if (fabs(matrix_computed(i, j) - matrix_expected(i, j)) > tolerance) {
                throw runtime_error("Sorry the elements of the matrices are not equal(with a certain tolerance), there is an error somewhere");
            }
        }
    }  
}
int main() {
    MatrixXd nodes=read_file("nodes.txt");
    MatrixXd triangles=read_file("elements.txt");
    MatrixXd boundaries=read_file("boundaries.txt");
    /*  3---4---5
        | / | / |
        2---9---6
        | / | / |
        1---8---7
        between each node there is a length of 1 to facilitat my calculation
        The following defined matrices are calculated thanks to the overleaf document with the matrices*/
    MatrixXd sigma(2, 2); // to facilite my calculations too
    sigma << 1, 0,
             0, 1;
    /*-----------Test for assembling matrix K------------*/
    MatrixXd computedK = compute_K(nodes, triangles, sigma);
    MatrixXd expected_k(9, 9);
    expected_k << 0.0005, -0.000167, 0, 0, 0, 0, 0, -0.00033, 0,
              -0.000167, 0.000833, -0.000167, 0, 0, 0, 0, 0, -0.0005,
              0, -0.000167, 0.00033, -0.000167, 0, 0, 0, 0, 0,
              0, 0, -0.000167, 0.001833, -0.00067, 0, 0, 0, -0.001,
              0, 0, 0, -0.00067, 0.0015, -0.000833333, 0, 0, 0,
              0, 0, 0, 0, -0.000833, 0.0033, -0.000833, 0, -0.00167,
              0, 0, 0, 0, 0, -0.000833333, 0.0015, -0.00067, 0,
              -0.00033, 0, 0, 0, 0, 0, -0.00067, 0.002, -0.001,
              0, -0.0005, 0, -0.001, 0, -0.00167, 0, -0.001, 0.004167;

    try {
        matrix_comparison(computedK, expected_k, 1e-4);// nb-ot to big tol because i didn't calculated it with a lot of terms after the coma
        cout << "Yihouuuuu the matrice K is correct" << endl;
    } catch (const runtime_error& e) {
        cerr << "Error: " << e.what() << endl;
    }
    /*-----------Test for assembling matrix f1------------*/
    MatrixXd computedf1= compute_f1(nodes,triangles,V_mu);
    MatrixXd expected_f1(9,9);

    expected_f1 << 5e-11, 8.33e-12, 0, 0, 0, 0, 0, 2.5e-11, 4.167e-11,
              8.33e-12, 6.667e-11, 8.333e-12, 4.1667e-11, 0, 0, 0, 0, 4.1667e-11,
              0, 8.333e-12, 1.667e-11, 1.667e-11, 0, 0, 0, 0, 0,
              0, 4.1667e-11, 1.6667e-11, 2.1667e-10, 5.833e-11, 0, 0, 0, 8.33e-11,
              0, 0, 0, 5.833e-11, 2.833e-10, 7.5e-11, 0, 0, 1.25e-10,
              0, 0, 0, 0, 7.5e-11, 3e-10, 7.5e-11, 0, 1.333e-10,
              0, 0, 0, 0, 0, 7.5e-11, 2.833e-10, 5.833e-11, 1.25e-10,
              2.5e-11, 0, 0, 0, 0, 0, 5.833e-11, 1.667e-10, 8.333e-11,
              4.1667e-11, 4.1667e-11, 0, 8.333e-11, 1.25e-10, 1.333e-10, 1.25e-10, 8.333e-11, 6.167e-10;
    try {
        matrix_comparison(computedf1, expected_f1, 1e-4);// nb-ot to big tol because i didn't calculated it with a lot of terms after the coma
        cout << "Yihouuuuu the matrice f1 is correct" << endl;
    } catch (const runtime_error& e) {
        cerr << "Error: " << e.what() << endl;
    }
    /*-----------Test for assembling matrix f2------------*/
    VectorXd computedf2=compute_f2(nodes,triangles,V_mfv);
    VectorXd expected_f2(9);
    expected_f2<< 8.333e-11,
              1.25e-10,
              4.167e-11,
              4.167e-10,
              5.4167e-10,
              5.833e-10,
              -8.333e-11,
              3.333e-10,
              8.75e-10;
    try {
        matrix_comparison(computedf2, expected_f2, 1e-4);// nb-ot to big tol because i didn't calculated it with a lot of terms after the coma
        cout << "Yihouuuuu the matrice f2 is correct" << endl;
    } catch (const runtime_error& e) {
        cerr << "Error: " << e.what() << endl;
    }
    /*-----------Test for assembling matrix H1------------*/
    MatrixXd computedH1_u=compute_H1(nodes,boundaries,rho_u);
    MatrixXd expected_H1u(9,9);
    expected_H1u<< 5.83333e-14,          0,          0,          0,          0,          0,          0, 5.83333e-14,          0,
                       0,          0,          0,          0,          0,          0,          0,          0,          0,
                       0,          0, 5.83333e-14, 5.83333e-14,          0,          0,          0,          0,          0,
                       0,          0, 5.83333e-14, 4.66667e-13,    1.75e-13,          0,          0,          0,          0,
                       0,          0,          0,    1.75e-13,    8.75e-13, 2.33333e-13,          0,          0,          0,
                       0,          0,          0,          0, 2.33333e-13, 9.33333e-13, 2.33333e-13,          0,          0,
                       0,          0,          0,          0,          0, 2.33333e-13,    8.75e-13,    1.75e-13,          0,
             5.83333e-14,          0,          0,          0,          0,          0,    1.75e-13, 4.66667e-13,          0,
                       0,          0,          0,          0,          0,          0,          0,          0,          0;
    try {
        matrix_comparison(computedH1_u, expected_H1u, 1e-4);// nb-ot to big tol because i didn't calculated it with a lot of terms after the coma
        cout << "Yihouuuuu the matrice H1_u is correct" << endl;
    } catch (const runtime_error& e) {
        cerr << "Error: " << e.what() << endl;
    }
    MatrixXd computedH1_v=compute_H1(nodes,boundaries,rho_v);
    MatrixXd expected_H1v(9,9);
    expected_H1v<<6.25e-14,        0,        0,        0,        0,        0,        0, 6.25e-14,        0,
                     0,        0,        0,        0,        0,        0,        0,        0,        0,
                     0,        0, 6.25e-14, 6.25e-14,        0,        0,        0,        0,        0,
                     0,        0, 6.25e-14,    5e-13, 1.875e-13,        0,        0,        0,        0,
                     0,        0,        0, 1.875e-13, 9.375e-13,  2.5e-13,        0,        0,        0,
                     0,        0,        0,        0,  2.5e-13,    1e-12,  2.5e-13,        0,        0,
                     0,        0,        0,        0,        0,  2.5e-13, 9.375e-13, 1.875e-13,        0,
             6.25e-14,        0,        0,        0,        0,        0, 1.875e-13,    5e-13,        0,
                     0,        0,        0,        0,        0,        0,        0,        0,        0;
    try {
        matrix_comparison(computedH1_v, expected_H1v, 1e-4);// nb-ot to big tol because i didn't calculated it with a lot of terms after the coma
        cout << "Yihouuuuu the matrice H1_v is correct" << endl;
    } catch (const runtime_error& e) {
        cerr << "Error: " << e.what() << endl;
    }
    /*-----------Test for assembling matrix H2------------*/
    VectorXd computedH2_u=compute_H2(nodes,boundaries,rho_v,Cu_amb);
    VectorXd expected_H2u =VectorXd::Zero(9);
    try {
        matrix_comparison(computedH2_u, expected_H2u, 1e-4);// nb-ot to big tol because i didn't calculated it with a lot of terms after the coma
        cout << "Yihouuuuu the matrice H2_u is correct" << endl;
    } catch (const runtime_error& e) {
        cerr << "Error: " << e.what() << endl;
    }
    MatrixXd computedH2_v=compute_H2(nodes,boundaries,rho_v,Cv_amb);
    VectorXd expected_H2v =VectorXd::Zero(9);
    try {
        matrix_comparison(computedH2_v, expected_H2v, 1e-4);// nb-ot to big tol because i didn't calculated it with a lot of terms after the coma
        cout << "Yihouuuuu the matrice H2_v is correct" << endl;
    } catch (const runtime_error& e) {
        cerr << "Error: " << e.what() << endl;
    }

    return 0;
}
