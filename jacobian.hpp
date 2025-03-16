#include <iostream>
#include "eigen-3.4.0/Eigen/Eigen"
#include "eigen-3.4.0/Eigen/Dense"
//#include "parameters.cpp"

using namespace std;
using namespace Eigen;
/*----Here there is a nonlinear version of the functions R and its second derivatives-----*/
namespace jacs{
VectorXd R_u(const VectorXd& Cu, const VectorXd& Cv) {
    return (V_mu * Cu.array())*((K_mu + Cu.array())*(1.+Cv.array()/K_mv)).inverse();
}
VectorXd R_v(const VectorXd& Cu, const VectorXd& Cv) {
    return r_q*R_u(Cu, Cv).array() +V_mfv*(1.+Cu.array()/K_mfu).inverse();
}
MatrixXd dRudu(const VectorXd& Cu, const VectorXd& Cv) {
    ArrayXd derivative = (V_mu*K_mu)*((1. + Cv.array()/K_mv)*(K_mu + Cu.array()).pow(2)).inverse();
    return derivative.matrix().asDiagonal();
}
MatrixXd dRudv(const VectorXd& Cu, const VectorXd& Cv) {
    ArrayXd derivative = (-K_mv * V_mu * Cu.array())*((K_mu + Cu.array()) * (K_mv + Cv.array()).pow(2)).inverse();
    return derivative.matrix().asDiagonal();
}
MatrixXd dRvdu(const VectorXd& Cu, const VectorXd& Cv) {
    ArrayXd derivative = r_q * (V_mu*K_mu)*((1. + Cv.array()/K_mv)*(K_mu + Cu.array()).pow(2)).inverse() - V_mfv*(K_mfu * (1. + Cu.array() / K_mfu).pow(2)).inverse();
    return derivative.matrix().asDiagonal();
}
MatrixXd dRvdv(const VectorXd& Cu, const VectorXd& Cv) {
    return r_q * dRudv(Cu, Cv);
}

}