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
/**
 * \file Tool.h
 * \brief File declaring convenient functions.

 * This tool file implements convenient functions for manipulating Eigen Containers,
 * import data from textfile,displaying containers, generate parameters and do some Math stuff.
 *
 */

#ifndef __Tool_H__
#define __Tool_H__
#include "Covfunc.h"


/**
 * \typedef TypePair
 * \brief typedef for the preference pairs.
 *
 * This typedef is used to store dynamics Eigen matrices of two cols into a dynamic Vector.
 */
typedef Matrix<Matrix<double, Dynamic, 2>, Dynamic, 1> TypePair;


/**
 * \typedef VecMat
 * \brief typedef for the matrices gradients covariances.
 *
 * This typedef is used to store Eigen dynamics matrices into a dynamic Vector.
 * Used into the project for storing the matrices gradients covariances.
 */
typedef Matrix<Matrix<double, Dynamic, Dynamic>, Dynamic, 1> VecMat;


/**
 * \typedef column_vector
 * \brief typedef for the dlib Vectors.
 *
 * This typedef is a dlib dynamic vector of double.
 */
typedef matrix<double, 0, 1> column_vector;

/**
 * \typedef  Vectbool
 * \brief Vector of boolean.
 *
 * This typedef is an Eigen dynamic vector of Boolean.
 */
//typedef of a vector containing boolean
typedef Matrix<bool, Dynamic, 1> Vectbool;




//input functions
    /*!
     *	\fn MatrixXd GetData(const string& myfile)
     *  \brief Gather a Matrix from a textfile into a MatrixXd container.
     *
     *	Gather a Matrix from a textfile. The first information into the textfile must be
     * the number of rows, then the number of column. The insertion into the Matrix is Row 
     * Major.
     *	\param string myfile
     */
MatrixXd GetData(const string& myfile);

    /*!
     *	\fn TypePair InputPair(const string& myfile)
     *  \brief Gather multiple Matrices from a textfile into a TypePair container.
     *
     *	Gather multiples Matrices from a folder. A text file names "Informations.txt"
     *	contains the numbers of matrices into the folder. Then this matrices are stored
     * into a TypePair container.
     *	\param string myfile
     */
TypePair InputPair(const string& myfile);

//convenient access functions

    /*!
     *	\fn fliplr(MatrixXd& a)
     *  \brief flip a Matrix.
     *
     *Inverse the parameters contained into the input parameter by row.
     *	\param MatrixXd a
     */
void fliplr(MatrixXd& a);
VectorXd find(const VectorXd a, double b);
void Add(VectorXd& a, double val);
MatrixXd GetMat(MatrixXd mat, VectorXd t1, VectorXd t2);
MatrixXd GetMatRow(MatrixXd mat, VectorXd t1);
void SetMatRow(MatrixXd& a, VectorXd& t1,MatrixXd &b);
VectorXd GetVec(VectorXd vec, VectorXd t1);
MatrixXd SetMatGenIdx(MatrixXd mat, VectorXd t1, VectorXd t2);
VectorXd GetMatGenIdx(MatrixXd mat, VectorXd t1);
int sub2ind(int dimrow, int dimcol, int row, int col);
VectorXd sub2ind(int dimrow, int dimcol, VectorXd setrow, VectorXd setcol);
VectorXd concatmat(const MatrixXd& a );
void GetTheta(VectorXd& theta_t, VectorXd& theta_x, double& sigma , VectorXd& theta, int dim_t, int dim_x);
VectorXd concatTheta(const VectorXd& theta_t, const  VectorXd& theta_x, double sigma);
column_vector EigentoDlib(VectorXd a);
VectorXd DlibtoEigen(column_vector a);

//Maths functions

double normcdf(double x);
VectorXd normcdf(VectorXd x);
double normpdf(double x);
VectorXd normpdf(VectorXd x);
MatrixXd Kron(MatrixXd mat1, MatrixXd mat2);
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
VectorXd concatsort(const VectorXd &a,const VectorXd& b);
MatrixXd make_query_toydata(TypePair Oracle, int query_idx, int test_idx);
void loss_query_toydata(double &loss, const MatrixXd& F, bool& stop, int test_user_idx, int best_item_idx);

//Debugging functions

void dsp(string s);
void dsp(MatrixXd a, string s);
void dsp(double a, string s);
void dspair(TypePair a, string txt);








#endif