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
#ifndef __CGppe_H__
#define __CGppe_H__
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <string>
#include "Covfunc.h"
using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;


class CGppe
{
public :
    CGppe(Covfunc *covfunc_t, Covfunc *covfunc_x);// Default Contructor
    CGppe(VectorXd fnew, MatrixXd Kxnew, MatrixXd Kinvnew, MatrixXd Wnew, MatrixXd Lnew); //Constructor with parameters
    CGppe(const CGppe & g); //Contructor by recopy
    ~CGppe(); // Destructor
    VectorXd Getf();
    MatrixXd GetKx();
    MatrixXd GetKinv();
    MatrixXd GetW();
    MatrixXd GetL();
    LLT<MatrixXd> Getllt();
    VectorXd Getmustar();
    VectorXd Getvarstar();
    double Getp();

    double maximum_expected_improvement(const VectorXd & theta_x, const VectorXd& theta_t, const double& sigma,
                                        const MatrixXd& t, const MatrixXd & x, const VectorXd& idx_global, const VectorXd& ind_t, const VectorXd& ind_x, MatrixXd tstar, int N, double fbest);

    double expected_voi(const VectorXd & theta_x, const VectorXd& theta_t, const double& sigma,
                        const MatrixXd& t, const MatrixXd & x, const TypePair& train_pairs, VectorXd& idx_global, VectorXd& ind_t, VectorXd& ind_x, MatrixXd test_pair, double fbest);


    void Elicit(const VectorXd & theta_x, const VectorXd& theta_t, const double& sigma, const MatrixXd& train_t, const MatrixXd &x, const TypePair & train_pairs
                , const MatrixXd & test_t, int test_user_idx, MatrixXd  idx_pairs, int  Maxiter, const TypePair & Oracle);

    void Make_Predictions_New_User(const VectorXd & theta_x, const VectorXd& theta_t, double& sigma, const MatrixXd& train_t, const MatrixXd &x, const TypePair & train_pairs,
                                   const VectorXd & idx_global, const VectorXd& idx_global_1, const VectorXd& idx_global_2,
                                   const VectorXd& ind_t, const VectorXd& ind_x, const MatrixXd & test_t, const MatrixXd& idx_pairs, const VectorXd& ftrue, const VectorXd& ytrue);

    void Predictive_Utility_Distribution(MatrixXd t, MatrixXd test_t, int N, VectorXd idx_global);

    void Predict_CGppe_Laplace(double sigma, MatrixXd t, MatrixXd x, VectorXd idx_global, VectorXd ind_t, VectorXd ind_x,
                              MatrixXd tstar, MatrixXd test_pair);

    void Approx_CGppe_Laplace(const VectorXd& theta_x, const VectorXd& theta_t,
                             const double& sigma, const MatrixXd& t, const MatrixXd& x, const TypePair &all_pairs, const VectorXd &idx_global, const VectorXd& idx_global_1, const VectorXd& idx_global_2, const VectorXd& ind_t, const VectorXd& ind_x, int M, int N);

    double log_likelihood(double sigma, TypePair all_pairs, VectorXd idx_global_1, VectorXd idx_global_2, int M, int N);

    VectorXd deriv_log_likelihood_CGppe_fast(double sigma, const TypePair & all_pairs, VectorXd idx_global_1, VectorXd idx_global_2, int M, int N);

    MatrixXd deriv2_log_likelihood_CGppe_fast(double sigma, TypePair all_pairs, VectorXd idx_global_1, VectorXd idx_global_2, int M, int N);

    double get_fbest(int N);
    //Variables
private :
    VectorXd f;
    MatrixXd Kx;
    MatrixXd Kinv;
    MatrixXd W;
    MatrixXd L;
    Covfunc *covfunc_x;
    Covfunc *covfunc_t;
    double p;
    VectorXd mustar;
    VectorXd varstar;
    LLT<MatrixXd> llt;
};
#endif