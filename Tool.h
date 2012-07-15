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

    /*!
     *	\fn VectorXd find(const VectorXd a, double b)
     *  \brief Return the indexes where is b.
     *
     *	This function returns a VectorXd which contains the indexes of the input
     *	VectorXd where the value b is stored.
     *	\param VectorXd a, double b
     */
VectorXd find(const VectorXd a, double b);

    /*!
     *	\fn VectorXd void Add(VectorXd& a, double val)
     *  \brief add a double to a VectorXd
     *
     *	This function returns add a row into the input VectorXd a and put inside
     *	the double val
     *	\param VectorXd a, double val
     */
void Add(VectorXd& a, double val);

    /*!
     *	\fn MatrixXd GetMat(MatrixXd mat, VectorXd t1, VectorXd t2)
     *  \brief Returns the selected elements of mat
     *
     *	Return a MatrixXd corresponding to the elements of mat given by
     *	the VectorXd t1 and t2.
     *	\param MatrixXd mat, VectorXd t1, VectorXd t2
     */
MatrixXd GetMat(MatrixXd mat, VectorXd t1, VectorXd t2);

    /*!
     *	\fn MatrixXd GetMatRow(MatrixXd mat, VectorXd t1)
     *  \brief Return the selected rows of mat.
     *
     *	Return a MatrixXd containing the rows of mat given by
     *	the VectorXd t1.
     *	\param MatrixXd mat, VectorXd t1
     */
MatrixXd GetMatRow(MatrixXd mat, VectorXd t1);

    /*!
     *	\fn void SetMatRow(MatrixXd& a, VectorXd& t1,MatrixXd &b)
     *  \brief Set the selected rows of a.
     *	Set the selected rows of the MatrixXd a given by t1 with the rows contained
     *	inside the MatrixXd b
     *	\param MatrixXd a, VectorXd t1, MatrixXd b
     */
void SetMatRow(MatrixXd& a, VectorXd& t1,MatrixXd &b);

    /*!
     *	\fn VectorXd GetVec(VectorXd vec, VectorXd t1)
     *  \brief Return the selected elements of vec.
     *
     *	Returns the selected elements of vec given by t1
     *	\param VectorXd vec, VectorXd t1
     */
VectorXd GetVec(VectorXd vec, VectorXd t1);

    /*!
     *	\fn MatrixXd SetMatGenIdx(MatrixXd mat, VectorXd t1, VectorXd t2)
     *  \brief Set the selected elements of mat with a general index.
     *
     *	Set the selected elements of mat with a general index given by t1
     * and the new value given by t2
     *	\param VectorXd mat, VectorXd t1, VectorXd t2
     */
MatrixXd SetMatGenIdx(MatrixXd mat, VectorXd t1, VectorXd t2);

    /*!
     *	\fn VectorXd GetMatGenIdx(MatrixXd mat, VectorXd t1)
     *  \brief Get the selected elements of mat with a general index.
     *
     *	Get the selected elements of mat with a general index given by t1
     *	\param VectorXd mat, VectorXd t1
     */
VectorXd GetMatGenIdx(MatrixXd mat, VectorXd t1);

    /*!
     *	\fn int sub2ind(int dimrow, int dimcol, int row, int col)
     *  \brief Convert substrict index to linear index.
     *
     *	Transform a substrict index into a linera index, by given
     *	the dimensions of te matrix and the 2 coordonates of the index.
     *	\param int dimrow, int dimcol, int row, int col
     */
int sub2ind(int dimrow, int dimcol, int row, int col);

    /*!
     *	\fn VectorXd sub2ind(int dimrow, int dimcol, VectorXd setrow, VectorXd setcol)
     *  \brief Convert a set of substrict index to linear index.
     *
     *	Transform a set of substrict index into a linera index, by given
     *	the dimensions of te matrix and the 2 coordonates of the index.
     *	\param int dimrow, int dimcol, VectorXd setrow, VectorXd setcol
     */
VectorXd sub2ind(int dimrow, int dimcol, VectorXd setrow, VectorXd setcol);

    /*!
     *	\fn VectorXd concatmat(const MatrixXd& a )
     *  \brief Convert a Matrix into a Vector
     *
     *	Fill a Vector with all the elements of te MatrixXd. The filling is column major.
     *	\param MatrixXd a
     */
VectorXd concatmat(const MatrixXd& a );

    /*!
     *	\fn void GetTheta(VectorXd& theta_t, VectorXd& theta_x, double& sigma , VectorXd& theta, int dim_t, int dim_x)
     *  \brief Set theta_x, theta_t and sigma
     *
     *	Set the hyperparameters theta_x, theta_t and sigma from a single Vector where
     *	all the informations where gathered. dim_t and dim_x give an indication about the parameters size.
     *	\param VectorXd theta_t, VectorXd theta_x, double sigma , VectorXd theta, int dim_t, int dim_x
     */
void GetTheta(VectorXd& theta_t, VectorXd& theta_x, double& sigma , VectorXd& theta, int dim_t, int dim_x);

    /*!
     *	\fn VectorXd concatTheta(const VectorXd& theta_t, const  VectorXd& theta_x, double sigma)
     *  \brief Concatenate the hyperparameters
     *
     *	Return a single Vector which contains theta_t, theta_x and sigma	
     *	\param VectorXd theta_t, VectorXd theta_x, double sigma , VectorXd theta, int dim_t, int dim_x
     */
VectorXd concatTheta(const VectorXd& theta_t, const  VectorXd& theta_x, double sigma);

    /*!
     *	\fn column_vector EigentoDlib(VectorXd a)
     *  \brief Convert an Eigen VectorXd into a dlib column_vector
     *
     *	\param VectorXd a
     */
column_vector EigentoDlib(VectorXd a);

    /*!
     *	\fn VectorXd DlibtoEigen(column_vector a)
     *  \brief Convert a dlib column_vector into a Eigen VectorXd
     *
     *	\param column_vector a
     */
VectorXd DlibtoEigen(column_vector a);

//Maths functions


    /*!
     *	\fn normcdf(double x)
     *  \brief Computes the normcdf
     *
     *	Computes the Normal cumulative distribution function with the given parameter x
     *	\param double x
     */
double normcdf(double x);

    /*!
     *	\fn VectorXd normcdf(VectorXd x)
     *  \brief Computes the normcdf from a set of parameters
     *
     *	Computes the Normal cumulative distribution function with each parameter into the VectorXd x and
     *	return a VectorXd containing the results.
     *	\param VectorXd x
     */
VectorXd normcdf(VectorXd x);

    /*!
     *	\fn normpdf(double x)
     *  \brief Computes the normpdf
     *
     *	Computes the Normal probability density function with the given parameter x
     *	\param double x
     */
double normpdf(double x);

    /*!
     *	\fn VectorXd normpdf(VectorXd x)
     *  \brief Computes the normcdf from a set of parameters
     *
     *	Computes the Normal probability density function  with each parameter into the VectorXd x and
     *	return a VectorXd containing the results.
     *	\param VectorXd x
     */
VectorXd normpdf(VectorXd x);

    /*!
     *	\fn MatrixXd Kron(MatrixXd mat1, MatrixXd mat2)
     *  \brief Computes the kronecker product
     *
     *	Returns a MatrixXd which is the result of the Kronecker product
     *	between mat1 and mat2.
     *	\param MatrixXd mat1, MatrixXd mat2
     */
MatrixXd Kron(MatrixXd mat1, MatrixXd mat2);

    /*!
     *	\fn VectorXd GetDiff(VectorXd a, VectorXd b)
     *  \brief Returns the difference between two Vectors
     *
     *	Returns a Vector filling with one or zero, regarding if
     *	The values contained into a and b are equals or not for the same index.
     *	\param VectorXd a, VectorXd b
     */
VectorXd GetDiff(VectorXd a, VectorXd b);

    /*!
     *	\fn VectorXd Get_Cumulative_Val(VectorXd idx, VectorXd val, int n)
     *  \brief Accumulates the values in val according to the indices given by idx
     *
     *	\param VectorXd idx, VectorXd val, int n
     */
VectorXd Get_Cumulative_Val(VectorXd idx, VectorXd val, int n);

    /*!
     *	\fn MatrixXd get_dWdf(VectorXd all_diag_idx, VectorXd f, int ind_t, int ind_x, double sigma, MatrixXd pairs, int  M, int N)
     *  \brief returns the Derivatives of the Matrix W wrt
     *
     *	Derivatives of the matrix W wrt a  single component in f. This is 
     *	useful in computing the implicit derivatives of the marginal
     *	likelihood (hyper-parameter learning).
     *	\param[in] f :  The current vector  f  (mode) 
     *	\param[in] ind_t :  scalar index of the current user
     *	\param[in] ind_x :  scalar index of the current item
     *	\param[in] sigma :  The noise parameter (scale factor of the preferences)
     *	\param[in] pairs :  Matrix of preference pairs for user indexed by ind_t 
     *	\param VectorXd all_diag_idx, VectorXd f, int ind_t, int ind_x, double sigma, MatrixXd pairs, int  M, int N
     */
MatrixXd get_dWdf(VectorXd all_diag_idx, VectorXd f, int ind_t, int ind_x, double sigma, MatrixXd pairs, int  M, int N);

    /*!
     *	\fn get_dsigma(MatrixXd& dWdsigma, double& dloglike_dsigma, const VectorXd& f, double sigma, const TypePair& all_pairs, int M, int N)
     *  \brief Returns the derivatives wrt sigma
     *
     *	Set the Matrix dWdsigma which is the Derivatives of the Matrix W wrt dtheta_sigma
     *	with theta_sigma = log(sigma). Set dloglike_dsigma which is the Derivative 
     *	of the conditional likelihood wrt theta_sigma
     *  . 
     *	\param[in] f :  The current vector  f  (mode) 
     *	\param[in] ind_t :  scalar index of the current user
     *	\param[in] ind_x :  scalar index of the current item
     *	\param[in] sigma :  The noise parameter (scale factor of the preferences)
     *	\param[in] all_pairs : TypePair of M elements. Each element is a O_m x 2 matrix 
     *  where O_m is the number of preferences observed for the corresponding
     *  user. Each row all_pairs(m) contains a preference relation 
     *  of the form all_pairs(m)(0) > all_pairs(m)(1) 
     *	\param[in] M : The number of users
     *	\param[in] N : The number of items  
     *	\param MatrixXd dWdsigma, double dloglike_dsigma, VectorXd f, double sigma, TypePair all_pairs, int M, int N
     */
void get_dsigma(MatrixXd& dWdsigma, double& dloglike_dsigma, const VectorXd& f, double sigma, const TypePair& all_pairs, int M, int N);

    /*!
     *	\fn VectorXd get_dlogp_dsigma(MatrixXd& dWdsigma, double& dloglike_dsigma, const VectorXd& f, double sigma, const TypePair& all_pairs, int M, int N)
     *  \brief Returns dfdsigma
     
     *	Returns the Derivative of the gradient of the conditional likelihood wrt sigma
     *	This is used to compute df_dsigma in the implicit derivative of
     *	the marginal likelihood wrt sigma. (hyper-parameter learning) 
     *  . 
     *	\param[in] f :  The current vector  f  (mode) 
     *	\param[in] ind_t :  scalar index of the current user
     *	\param[in] ind_x :  scalar index of the current item
     *	\param[in] sigma :  The noise parameter (scale factor of the preferences)
     *	\param[in] all_pairs : TypePair of M elements. Each element is a O_m x 2 matrix 
     *  where O_m is the number of preferences observed for the corresponding
     *  user. Each row all_pairs(m) contains a preference relation 
     *  of the form all_pairs(m)(0) > all_pairs(m)(1) 
     *	\param[in] M : The number of users
     *	\param[in] N : The number of items 
     *	\param[out] N : dfdsigma: vector of derivatives d_dsigma(d_df (log p (D | f) ) )

     *	\param MatrixXd dWdsigma, double dloglike_dsigma, VectorXd f, double sigma, TypePair all_pairs, int M, int N
     */
VectorXd get_dlogp_dsigma(MatrixXd& dWdsigma, double& dloglike_dsigma, const VectorXd& f, double sigma, const TypePair& all_pairs, int M, int N);

    /*!
     *	\fn unsigned long fact(int num)
     *  \brief Returns the Factorial of the input number
     *
     *	\param int num
     */
unsigned long fact(int num);

    /*!
     *	\fn BinCoef(int n, int k)
     *  \brief Returns Binomial coefficient C(n, k)
     *
     *	\param int n, int 
     */
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