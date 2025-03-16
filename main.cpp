#include <iostream>
#include "eigen-3.4.0/Eigen/Eigen"
#include "eigen-3.4.0/Eigen/Dense"
#include "matrix.cpp"
#include "newton_raphson.cpp"
//#include "parameters.cpp"
#include "read_file.hpp"
#include "jacobian.hpp"
#include <fstream>
#include <chrono>
#include <string>
#include <stdio.h>
#include <stdlib.h>

using namespace std;
using namespace Eigen;


/*Run the code with *  g++ -o test main.cpp
                    * ./test <number between 1 and 6> */

int main(int argc, char *argv[]){
    
  /*----Load constants from parameters.cpp----*/
  // Parameters depend on storage conditions:
  // storage = 1 if Orchard
  // storage = 2 if Shelf life
  // storage = 3 if Refrigerator
  // storage = 4 if Precooling
  // storage = 5 if Disorder inducing
  // storage = 6 if Optimal CA
  /*-----We will ask the user which storage conditions he want------*/
  int storage;

  if (argc!=2) { 
      printf("Launch it like this : ./test <storage_number>\n");
      printf("With a choosen <storage_number>: \n");
      printf("1: Orchard, 2: Shelf life, 3: Refrigerator, 4: Precooling, 5: Disorder inducing, 6: Optimal CA\n");
      return 1;
  }

  storage=atoi(argv[1]);

  if (storage<1 || storage>6) {
      printf("Please choose a storage number between 1 and 6.\n");
      printf("Remember -> 1: Orchard, 2: Shelf life, 3: Refrigerator, 4: Precooling, 5: Disorder inducing, 6: Optimal CA\n");
      return 1; 
  }

  auto start = chrono::high_resolution_clock::now();
  setParameters(storage);
  

  /*-------Read files to store informations---------*/

  MatrixXd nodes= read_file("../semi_circle/readNodes.txt");
  MatrixXd triangles= read_file("../semi_circle/readElements.txt");
  MatrixXd boundaries= read_file("../semi_circle/readBoundaries.txt");
  

  /*----------Assemble the system from matrix.cpp------------*/
  
  /*---First part of the system-----*/
  MatrixXd sigma_u = Matrix<double, 2, 2>::Zero();
  MatrixXd sigma_v = Matrix<double, 2, 2>::Zero();
  sigma_u(0,0) = sig_ur;
  sigma_u(1,1) = sig_uz;
  sigma_v(0,0) = sig_vr;
  sigma_v(1,1) = sig_vz;

  MatrixXd Ku=compute_K(nodes,triangles,sigma_u);
  MatrixXd Kv=compute_K(nodes,triangles,sigma_v);
  if(Ku.rows()!= nodes.rows() || Ku.cols() != nodes.rows() || Kv.rows() != nodes.rows() || Kv.cols() != nodes.rows()) {
          throw std::runtime_error("Matrix Ku or Kv of wrong dimension");
  }
  ofstream file_j;
  file_j.open("results/result_Kv.txt");
  file_j << Kv;
  file_j.close();

  /*----Second part of the system-----*/

  MatrixXd fu = -compute_f1(nodes,triangles,V_mu);
  MatrixXd fv1 = compute_f1(nodes,triangles,V_mu);
  MatrixXd fv2 = compute_f2(nodes,triangles,V_mfv);
  if(fu.rows()!= nodes.rows() || fu.cols() != nodes.rows() || fv1.rows() != nodes.rows() || fv1.cols() != nodes.rows()) {
          throw std::runtime_error("Matrix fu or fv1 of wrong dimension");
  }
  if(fv2.rows()!= nodes.rows() || fv2.cols() != 1) {
          throw std::runtime_error("Matrix fu or fv1 of wrong dimension");
  }
  /*-----Third part of the system------*/

  MatrixXd Hu1 = compute_H1(nodes,boundaries, rho_u);
  MatrixXd Hu2 = compute_H2(nodes,boundaries, rho_u, Cu_amb);
  MatrixXd Hv1 = compute_H1(nodes,boundaries, rho_v);
  MatrixXd Hv2 = compute_H2(nodes,boundaries, rho_v, Cv_amb);
  if(Hu1.rows()!= nodes.rows() || Hu1.cols() != nodes.rows() || Hv1.rows() != nodes.rows() || Hv1.cols() != nodes.rows()) {
          throw std::runtime_error("Matrix Hu1 or Hv1 of wrong dimension");
  }
  if(Hu2.rows()!= nodes.rows() || Hu2.cols() != 1 || Hv2.rows() != nodes.rows() || Hv2.cols() != 1) {
          throw std::runtime_error("Matrix Hu2 or Hv2 of wrong dimension");
  }

  /*------Assemble the final system-----*/

  int M = nodes.rows();
  VectorXd c_u = VectorXd::Zero(M);
  VectorXd c_v = VectorXd::Zero(M);

//   fu_linear = fu*V_mu/K_mu*c_u;
//   fu_nonlin = fu*R_u(c_u, c_v);

//   fv_linear = fv1*r_q*c_u + V_mfv*fv2;
//   fv_nonlin = fv1*R_v(c_u, c_v);

  // Combined function for both equations
  auto F_combined = [&] (const VectorXd& x) {
    int M = Ku.rows(); 
    VectorXd c_u = x.head(M);
    VectorXd c_v = x.tail(M);
    VectorXd eq1 = Ku*c_u - fu*jacs::R_u(c_u, c_v) + Hu1*c_u - Hu2;
    VectorXd eq2 = Kv*c_v - fv1*jacs::R_v(c_u, c_v) + Hv1*c_v - Hv2;
    VectorXd combined = VectorXd(2 * M);
    combined << eq1, eq2;
    return combined;
  };
  
  /*----------Computation of the Jacobian-----------*/
  auto J_combined = [&] (const VectorXd& x) {
    int M = x.rows() ; 
    VectorXd c_u = x.head(M/2);
    VectorXd c_v = x.tail(M/2);

    MatrixXd dR_udu = jacs::dRudu(c_u, c_v);
    MatrixXd dR_udv = jacs::dRudv(c_u, c_v);
    MatrixXd dR_vdu = jacs::dRvdu(c_u, c_v);
    MatrixXd dR_vdv = jacs::dRvdv(c_u, c_v); 

    MatrixXd J(M, M);
    J << Ku - fu*dR_udu + Hu1, -fu*dR_udv, -fv1*dR_vdu, Kv - fv1*dR_vdv + Hv1; // Derivation of the system equations in c_u and c_v

    return J;
  };  
  

  /*------Using solver i.e. Newton Raphson------*/

  VectorXd x0 = VectorXd::Zero(2*nodes.rows());
  
  // To initialize x0, solve the linear system.
  // Ku*c_u - fu*V_mu/K_mu*c_u + Hu1*c_u - Hu2 = 0
  // (Ku - fu*V_mu/K_mu + Hu1)*c_u = Hu2
  VectorXd c_u0 = (Ku - fu*V_mu/K_mu + Hu1).colPivHouseholderQr().solve(Hu2);
  
  // Kv*c_v - fv1*r_q*c_u - V_mfv*fv2 + Hv1*c_v - Hv2 = 0
  // (Kv + Hv1)*c_v = fv1*r_q*c_u0 + V_mfv*fv2 + Hv2
  VectorXd c_v0 = (Kv + Hv1).colPivHouseholderQr().solve(fv1*r_q*c_u0 + V_mfv*fv2 + Hv2);

  x0 << c_u0, c_v0;

  /*-------Verify the correctness of our jacobian i.e. using finite differences------*/
  // we will choose a specific point i.e. x0 and apply finite difference formula with a certain step
  // M=Ku.rows();
  // MatrixXd J_approx = MatrixXd::Zero(2*M, 2*M);
  // double h=1e-5; 
  // VectorXd x_plus_h;
  // VectorXd f_x=F_combined(x0); 
  // for(int j =0; j <2*M; ++j) {
  //     x_plus_h=x0;
  //     x_plus_h(j)+=h; 
  //     VectorXd f_x_plus_h=F_combined(x_plus_h);
  //     J_approx.col(j)=(f_x_plus_h-f_x)/h;
  // }
  // MatrixXd J = J_combined(x0);
  // double norm_diff = (J - J_approx).norm();
  // cout << "The difference in norm between J and J_approx at point x0 is: " << norm_diff << endl;
  
  
  ofstream file;
  file.open("results/result_c0.txt");
  file << x0;
  file.close();
  
  newton_raphson(F_combined, J_combined, x0);
  if(x0.rows()!= 2*nodes.rows() || x0.cols() != 1 ) {
            throw std::runtime_error("Vector x0 after Newton Raphson not the right dimension anymore");
    }

  /* -------- Write results to files --------------- */
  M = x0.rows() ; 
  c_u = x0.head(M/2);
  c_v = x0.tail(M/2);
  
  ofstream file_u;
  file_u.open("results/result_cu.txt");
  file_u << c_u;
  file_u.close();

  ofstream file_v;
  file_v.open("results/result_cv.txt");
  file_v << c_v;
  file_v.close();
  
  /*---------Timing to compare performances with matlab-----------*/
  auto stop = chrono::high_resolution_clock::now(); 
  auto duration = chrono::duration_cast<chrono::seconds>(stop-start);
  std::cout << "It took "<< duration.count() <<" seconds to find a solution reaching the defined tolerance.\n";
  return 0;


}
