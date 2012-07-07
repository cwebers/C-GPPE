// Copyright (c) 2012, National ICT Australia
// All rights reserved.
//
// The contents of this file are subject to the Mozilla Public License
// Version 2.0 (the "License"); you may not use this file except in
// compliance with the License. You may obtain a copy of the License at
// http://www.mozilla.org/MPL/
//
// Software distributed under the License is distributed on an "AS IS"
// basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
// License for the specific language governing rights and limitations
// under the License.
#ifndef __Learn_H__
#define __Learn_H__
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <string>
#include "Gppe.h"
using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;

class CLearner
{
public :
    CLearner();
    
    CLearner(Covfunc *Covt, Covfunc *Covx,
          MatrixXd T, MatrixXd X, TypePair All_pairs, VectorXd Idx_global, VectorXd Idx_global_1,
          VectorXd Idx_global_2, VectorXd Ind_t, VectorXd Ind_x, int  m, int n);// Default Contructor
    //Learn(const Learn & l); //Contructor by recopy

    ~CLearner(); // Destructor


    //  const double negative_marginal_log_likelihood(const column_vector &dltheta);
    const double negative_marginal_log_likelihood(const column_vector & dltheta);
    
    //  VectorXd gradient_negative_marginal_loglikelihood(VectorXd theta);
    column_vector gradient_negative_marginal_loglikelihood(const column_vector & theta);

    //  Calculate function at point
    // double operator() ( const column_vector & point) const;

    //  Calculate gradient at point
    // column_vector & operator() ( const column_vector& point) const;
    // int operator() ( const column_vector & point) const;

// {
// VectorXd theta=DlibtoEigen(arg);
// VectorXd theta_x, theta_t;
// double sigma;
// GetTheta(theta_x, theta_t, sigma, theta);
// sigma=exp(sigma);
// covt->SetTheta(theta_t);
// covx->SetTheta(theta_x);
// Gppe g = Gppe(covt, covx);
// g.Approx_Gppe_Laplace( theta_x, theta_t, sigma,
//    t, x, train_pairs, idx_global, idx_global_1, idx_global_2, ind_t, ind_x, M, N);
//
//    double cond_loglike=g.log_likelihood(sigma, train_pairs, idx_global_1, idx_global_2, M, N);
// VectorXd fvis=GetVec(g.Getf(),idx_global);
// double margl=-0.5*(-log(g.GetKinv().determinant())+2*(log(g.GetL().diagonal().array()).sum()))
// -0.5*fvis.transpose()*g.GetKinv()*fvis +cond_loglike;
// dsp(-margl,"nl");
// return -margl;
//        return this->negative_marginal_log_likelihood(arg);
//}
    //VectorXd Optimize(VectorXd theta_first);
    //Variables
private :
    Covfunc *covx;
    Covfunc *covt;
    MatrixXd t;
    MatrixXd x;
    VectorXd idx_global;
    VectorXd idx_global_1;
    VectorXd idx_global_2;
    VectorXd ind_t;
    VectorXd ind_x;
    TypePair train_pairs;
    int M;
    int N;


};
#endif