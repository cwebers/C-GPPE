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
#include "Covfunc.h"



// typdef for the preference pairs : dynamic vectors of dynamic matrices
typedef Matrix<Matrix<double, Dynamic, 2>, Dynamic, 1> TypePair;

//typedef used for storing the matrices gradients covariances
typedef Matrix<Matrix<double, Dynamic, Dynamic>, Dynamic, 1> VecMat;

//typedef for the vector type in dlib
typedef matrix<double, 0, 1> column_vector;

//typedef of a vector containing boolean
typedef Matrix<bool, Dynamic, 1> Vectbool;



//no classed yet
MatrixXd make_query_toydata(TypePair Oracle, int query_idx, int test_idx);
void loss_query_toydata(double &loss, const MatrixXd& F, bool& stop, int test_user_idx, int best_item_idx);

//input functions
int GetDataline(const string& myfile);
int GetDatacol(const string& myfile);
MatrixXd GetData(const string& myfile);
TypePair InputPair(const string& myfile);

//convenient access functions

void fliplr(MatrixXd& a);
//int find(const MatrixXd& a, double val );
VectorXd find(const VectorXd a, int b);
void Add(VectorXd& a, double val);
//VectorXd find(const VectorXd& a, const VectorXd& b);
MatrixXd GetMat(MatrixXd mat, VectorXd t1, VectorXd t2);
MatrixXd GetMatRow(MatrixXd mat, VectorXd t1);
void SetMatRow(MatrixXd& a, VectorXd& t1,MatrixXd &b);
VectorXd GetVec(VectorXd vec, VectorXd t1);
MatrixXd SetMatGenIdx(MatrixXd mat, VectorXd t1, VectorXd t2);
VectorXd GetMatGenIdx(MatrixXd mat, VectorXd t1);
int sub2ind(int dimrow, int dimcol, int row, int col);
VectorXd sub2ind(int dimrow, int dimcol, VectorXd setrow, VectorXd setcol);
VectorXd concatmat(const MatrixXd& a );
void GetTheta(VectorXd& theta_x, VectorXd& theta_t, double& sigma , VectorXd& theta, int dim_t, int dim_x);
VectorXd concatTheta(const VectorXd& theta_t, const  VectorXd& theta_x, double sigma);
column_vector EigentoDlib(VectorXd a);
VectorXd DlibtoEigen(column_vector a);

//Maths functions

double normcdf(double x);
VectorXd normcdf(VectorXd x);
double normpdf(double x);
VectorXd normpdf(VectorXd x);
MatrixXd Kron(MatrixXd mat1, MatrixXd mat2);
MatrixXd SetNaN(MatrixXd a);
VectorXd MyNaNMean(MatrixXd a);
VectorXd GetDiff(VectorXd a, VectorXd b);
VectorXd Get_Cumulative_Val(VectorXd idx, VectorXd val, int n);
MatrixXd get_dWdf(VectorXd all_diag_idx, VectorXd f, int ind_t, int ind_x, double sigma, MatrixXd pairs, int  M, int N);
VectorXd get_cum2(VectorXd idx, VectorXd val, int n);
void get_dsigma(MatrixXd& dWdsigma, double& dloglike_dsigma, const VectorXd& f, double sigma, const TypePair& all_pairs, int M, int N);
VectorXd get_dlogp_dsigma(MatrixXd& dWdsigma, double& dloglike_dsigma, const VectorXd& f, double sigma, const TypePair& all_pairs, int M, int N);
unsigned long fact(int num);
unsigned long BinCoef(int n, int k);



//Generating parameters functions
void Generate(MatrixXd &idx_pairs,MatrixXd & t,MatrixXd& x,TypePair & Oracle, TypePair & train_pairs, MatrixXd & F, Covfunc *covx, Covfunc *covt, VectorXd &theta_t, VectorXd& theta_x, int &M, int &N,
VectorXd &ftrue, VectorXd &ytrue, MatrixXd& test_pairs, MatrixXd& test_t, MatrixXd &train_t );
void compute_global_index(VectorXd& idx_global_1, VectorXd& idx_global_2, const TypePair& all_pairs, int N);
void unique(VectorXd& a, const VectorXd& b, const VectorXd& c);
void ind2sub(VectorXd& ind_i, VectorXd& ind_j, int dimrow, int dimcol, VectorXd idx );
VectorXd ind2global(VectorXd vec, int j, int N);
VectorXd ind2global(VectorXd a, VectorXd b, int N);
VectorXd Nfirst(int N);
MatrixXd reshape(VectorXd f, int a, int b);
MatrixXd MatAdd(MatrixXd mat, MatrixXd ln);
MatrixXd combnk(int n);
VectorXd randperm(int n);


//Debugging functions

void dsp(string s);
void dsp(MatrixXd a, string s);
void dsp(double a, string s);
void dspair(TypePair a, string txt);








#endif