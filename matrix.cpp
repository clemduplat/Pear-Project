#include <iostream>
#include "eigen-3.4.0/Eigen/Eigen"
#include "eigen-3.4.0/Eigen/Dense"
#include "parameters.cpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;

//----------First part of the system----------//
MatrixXd compute_K(MatrixXd& nodes, MatrixXd& triangles,MatrixXd& sigma){
    MatrixXd K= MatrixXd::Zero(nodes.rows(),nodes.rows()); 
    Eigen::Matrix<double,3,2> coordinate;
    Eigen::Matrix<double,3,3> contribution;
    //cout << triangles.rows();
    for(unsigned i=0;i<triangles.rows();i++){
        int pt_1= triangles(i,0)-1;
        int pt_2= triangles(i,1)-1;
        int pt_3= triangles(i,2)-1;
        if((pt_1 < 0) || (pt_1 >= nodes.rows()) || (pt_2 < 0) || (pt_2 >= nodes.rows()) || (pt_3 < 0) || (pt_3 >= nodes.rows())) {
            std::cerr << "Index out of bounds at : " << pt_1 << ", " << pt_2 << ", " << pt_3 << std::endl;
            throw std::out_of_range("Node index out of bounds");
        }

        double area=triangles(i,3);
        // coord (r,z)
        double r_1=nodes(pt_1,0);
        double r_2=nodes(pt_2,0);
        double r_3=nodes(pt_3,0);
        double z_1=nodes(pt_1,1);
        double z_2=nodes(pt_2,1);
        double z_3=nodes(pt_3,1);
        //first matrix
        coordinate(0,0)=z_2-z_3;
        coordinate(1,0)=z_3-z_1;
        coordinate(2,0)=z_1-z_2;
        coordinate(0,1)=r_3-r_2;
        coordinate(1,1)=r_1-r_3;
        coordinate(2,1)=r_2-r_1;
        // matrix contribution of one element
        if(pt_1 < 0 || pt_1 >= K.rows() || pt_2 < 0 || pt_2 >= K.rows() || pt_3 < 0 || pt_3 >= K.rows()) {
            std::cerr << "K matrix index out of bounds." << std::endl;
            throw std::out_of_range("Node index out of bounds for K matrix");
        }

        contribution= ((r_1+r_2+r_3)/(12.*area))*coordinate*sigma*coordinate.transpose();
        //add the contribution to the matrix K for each element
        K(pt_1, pt_1) += contribution(0, 0);
        K(pt_2, pt_1) += contribution(1, 0);
        K(pt_1, pt_2) += contribution(1, 0);
        K(pt_3, pt_1) += contribution(2, 0);
        K(pt_1, pt_3) += contribution(2, 0);
        K(pt_2, pt_2) += contribution(1, 1);
        K(pt_2, pt_3) += contribution(1, 2);
        K(pt_3, pt_2) += contribution(1, 2);
        K(pt_3, pt_3) += contribution(2, 2);
        
    
    }
    return K;
}

//----------Second part of the system----------//
MatrixXd compute_f1(MatrixXd& nodes, MatrixXd& triangles, double V_mu){
    MatrixXd f1= MatrixXd::Zero(nodes.rows(),nodes.rows()); 
    Eigen::Matrix<double,3,3> contribution;
    for(unsigned i=0;i<triangles.rows();i++){
        int pt_1= triangles(i,0)-1;
        int pt_2= triangles(i,1)-1;
        int pt_3= triangles(i,2)-1;
        double area=triangles(i,3);

        if((pt_1 < 0) || (pt_1 >= nodes.rows()) || (pt_2 < 0) || (pt_2 >= nodes.rows()) || (pt_3 < 0) || (pt_3 >= nodes.rows())) {
            std::cerr << "Index out of bounds at : " << pt_1 << ", " << pt_2 << ", " << pt_3 << std::endl;
            throw std::out_of_range("Node index out of bounds");
        }
        
        double r_1=nodes(pt_1,0);
        double r_2=nodes(pt_2,0);
        double r_3=nodes(pt_3,0);

        contribution(0, 0)=6.*r_1+2.*r_2+2.*r_3;
        contribution(0, 1)=2.*r_1+2.*r_2+r_3;
        contribution(0, 2)=2.*r_1+r_2+2.*r_3;
        contribution(1, 0)=contribution(0, 1);
        contribution(1, 1)=2.*r_1+6.*r_2+2.*r_3;
        contribution(1, 2)=r_1+2.*r_2+2.*r_3;
        contribution(2, 0)=contribution(0, 2);
        contribution(2, 1)=contribution(1, 2);
        contribution(2, 2)=2.*r_1+2.*r_2+6.*r_3;

        // matrix contribution of one element
        contribution *= area/60.;

        //add the contribution to the matrix f1 for each element
        f1(pt_1, pt_1) += contribution(0, 0);
        f1(pt_1, pt_2) += contribution(0, 1);
        f1(pt_1, pt_3) += contribution(0, 2);
        f1(pt_2, pt_1) += contribution(1, 0);
        f1(pt_2, pt_2) += contribution(1, 1);
        f1(pt_2, pt_3) += contribution(1, 2);
        f1(pt_3, pt_1) += contribution(2, 0);
        f1(pt_3, pt_2) += contribution(2, 1);
        f1(pt_3, pt_3) += contribution(2, 2);
    
    }
    return f1;

}

VectorXd compute_f2(MatrixXd& nodes, MatrixXd& triangles, double V_mfv){
    VectorXd f2= VectorXd::Zero(nodes.rows()); 
    Eigen::Matrix<double,3,1> contribution;
    for(unsigned i=0;i<triangles.rows();i++){
        int pt_1= triangles(i,0)-1;
        int pt_2= triangles(i,1)-1;
        int pt_3= triangles(i,2)-1;
        double area=triangles(i,3);

        if((pt_1 < 0) || (pt_1 >= nodes.rows()) || (pt_2 < 0) || (pt_2 >= nodes.rows()) || (pt_3 < 0) || (pt_3 >= nodes.rows())) {
            std::cerr << "Index out of bounds at : " << pt_1 << ", " << pt_2 << ", " << pt_3 << std::endl;
            throw std::out_of_range("Node index out of bounds");
        }

        double r_1=nodes(pt_1,0);
        double r_2=nodes(pt_2,0);
        double r_3=nodes(pt_3,0);

        contribution(0, 0)=r_3-r_1;
        contribution(1, 0)=r_1+2.*r_2+r_3;
        contribution(2, 0)=r_1+r_2+2.*r_3;

        contribution *= area/12.;

        f2(pt_1) += contribution(0, 0);
        f2(pt_2) += contribution(1, 0);
        f2(pt_3) += contribution(2, 0);

    }
    return f2;

}

//----------Third part of the system----------//
MatrixXd compute_H1(MatrixXd& nodes, MatrixXd& boundaries, double rho){
    MatrixXd H1= MatrixXd::Zero(nodes.rows(),nodes.rows()); 
    Eigen::Matrix<double,2,2> contribution;
    for(unsigned i=0;i<boundaries.rows();i++){
        int pt_1= boundaries(i,0)-1;
        int pt_2= boundaries(i,1)-1;

        if((pt_1 < 0) || (pt_1 >= nodes.rows()) || (pt_2 < 0) || (pt_2 >= nodes.rows())) {
            std::cerr << "Index out of bounds at : " << pt_1 << ", " << pt_2 << std::endl;
            throw std::out_of_range("Node index out of bounds");
        }

        double r_1=nodes(pt_1,0);
        double r_2=nodes(pt_2,0);
        double z_1=nodes(pt_1,1);
        double z_2=nodes(pt_2,1);

        double length=boundaries(i, 2);

        contribution(0, 0)=3.*r_1+r_2;
        contribution(0, 1)=r_1+r_2;
        contribution(1, 0)=contribution(0, 1);
        contribution(1, 1)=r_1+3.*r_2;

        contribution *= length;

        H1(pt_1, pt_1) += contribution(0, 0);
        H1(pt_1, pt_2) += contribution(0, 1);
        H1(pt_2, pt_1) += contribution(1, 0);
        H1(pt_2, pt_2) += contribution(1, 1);
    
    }

    H1 *= (rho/12.);
    return H1;
}


VectorXd compute_H2(MatrixXd& nodes, MatrixXd& boundaries, double rho, double C_amb){
    VectorXd H2= VectorXd::Zero(nodes.rows()); 
    Eigen::Matrix<double,2,1> contribution;
    for(unsigned i=0;i<boundaries.rows();i++){
        int pt_1= boundaries(i,0)-1;
        int pt_2= boundaries(i,1)-1;

        if((pt_1 < 0) || (pt_1 >= nodes.rows()) || (pt_2 < 0) || (pt_2 >= nodes.rows())) {
            std::cerr << "Index out of bounds at : " << pt_1 << ", " << pt_2 << std::endl;
            throw std::out_of_range("Node index out of bounds");
        }

        double r_1=nodes(pt_1,0);
        double r_2=nodes(pt_2,0);
        double z_1=nodes(pt_1,1);
        double z_2=nodes(pt_2,1);

        double length=boundaries(i, 2);

        contribution(0, 0)=2.*r_1+r_2;
        contribution(1, 0)=r_1+2.*r_2;

        contribution *= length;

        H2(pt_1) += contribution(0, 0);
        H2(pt_2) += contribution(1, 0);
    
    }

    H2 *= (rho*C_amb/6.);
    return H2;
}

