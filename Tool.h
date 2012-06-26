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
#ifndef __Tool_H__
#define __Tool_H__
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <string>
#include <cmath>
using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::Dynamic;
using Eigen::Matrix;

typedef Matrix<Matrix<double,Dynamic,2>,Dynamic,1> TypePair;
VectorXd ind2global(VectorXd vec,int j,int N);
MatrixXd GetMat(MatrixXd mat,VectorXd t1, VectorXd t2);
VectorXd GetVec(VectorXd vec,VectorXd t1);
double normcdf(double x);
VectorXd normcdf(VectorXd x);
double normpdf(double x);
VectorXd normpdf(VectorXd x);
VectorXd Get_Cumulative_Val(VectorXd idx, VectorXd val, int n);
int sub2ind(VectorXd dim, int row, int col);
VectorXd sub2ind(int dimrow,int dimcol, VectorXd setrow,VectorXd setcol);
MatrixXd SetMatGenIdx(MatrixXd mat,VectorXd t1, VectorXd t2);
VectorXd GetMatGenIdx(MatrixXd mat,VectorXd t1);



#endif