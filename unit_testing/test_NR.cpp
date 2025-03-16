#include <iostream>
#include "./eigen-3.4.0/Eigen/Eigen"
#include "./eigen-3.4.0/Eigen/Dense"
#include "./newton_raphson.cpp"  

using Eigen::VectorXd;
/* This is a test function for the newton Raphson solver*/
/*verified using https://www.codesansar.com/numerical-methods/newton-raphson-method-online-calculator.htm */
//function
VectorXd F(const VectorXd& x) {
    VectorXd result(1);
    //result[0] = x[0] * x[0] - 4; 
    result[0]= 3*x[0] - cos(x[0]) - 1;
    return result;
}

//jacobian -> derivative
VectorXd J(const VectorXd& x) {
    VectorXd result(1);
    //result[0] = 2 * x[0]; 
    result[0]=3 + sin(x[0]);
    return result;
}

int main() {
    VectorXd x0(1);
    //Initial guess
    x0[0] = 2.0; 
    // by default 1e-6
    double tol = 1e-6; //0.000001
    newton_raphson(F, J, x0, tol);
    std::cout << "Root: " << x0[0] << std::endl;
    return 0;
}
