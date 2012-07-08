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
#include "CGppe.h"
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
    
      //VectorXd gradient_negative_marginal_loglikelihood(VectorXd theta);
    column_vector gradient_negative_marginal_loglikelihood(const column_vector & point);


    
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