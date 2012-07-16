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
#include "Tool.h"

/*! \class CGppe
   * \brief Main module used to compute Prediction and Elicitation
   *
   *	This Class Computes all the essentials methods for the Prediction and Elicitation processes.
   *	The members it owns are redundant variables needed for many functions into this module.
   */
class CGppe
{
public :

   /*! 
   * \brief Constructor which instanciates the covariances used for the process
   *
   */
    CGppe(Covfunc *covfunc_t, Covfunc *covfunc_x);
    
   /*! 
   * \brief Constructor with instanciation of the members
   *
   */
    CGppe(VectorXd fnew, MatrixXd Kxnew, MatrixXd Kinvnew, MatrixXd Wnew, MatrixXd Lnew); //Constructor with parameters
    
   /*! 
   * \brief Constructor by recopy
   *
   */
    CGppe(const CGppe & g);
    
   /*! 
   * \brief Destructor
   *
   */
    ~CGppe();
    
   /*! 
   * \brief Accessor for f
   *
   */
    VectorXd Getf();
    
   /*! 
   * \brief Accessor for Kx
   *
   */
    MatrixXd GetKx();
    
   /*! 
   * \brief Accessor for Kinv
   *
   */
    MatrixXd GetKinv();
    
   /*! 
   * \brief Accessor for W
   *
   */
    MatrixXd GetW();
    

   /*! 
   * \brief Accessor for L
   *
   */
    MatrixXd GetL();

   /*! 
   * \brief Accessor for llt
   *
   */
    LLT<MatrixXd> Getllt();
    
   /*! 
   * \brief Accessor for mustar
   *
   */
    VectorXd Getmustar();
    
   /*! 
   * \brief Accessor for varstar
   *
   */
    VectorXd Getvarstar();
    
   /*! 
   * \brief Accessor for p
   *
   */
    double Getp();
    
    
   /*! 
   * \brief Computes the [mei]
   *
   *	Computes the maximum expected improvement (MEI) of recommending items
   *	given by indices idx_xstar on user with features tstar.
   *
   *	\param[in] ind_t :  Indices of seen tasks.
   *	\param[in] ind_x :  Indices of seen items.
   *	\param[in] theta_x : The parameters for the items's covariance.
   *	\param[in] theta_t : The parameters for the users's covariance.
   *	\param[in] sigma : The noise parameter.
   *	\param[in] t :  The Matrix of User's features. 
   *	\param[in] x :  he Matrix of Item's features.
   *	\param[in] train_pairs : The Matrices of preferences pairs containing all the the preferences pairs of the users
   *	excepted the preferences pairs of the tested user.
   *	\param[in] idx_global : The unique global indices of the observed preferences.
   *	\param[in] tstar : Features of test user. 
   *	\param[in] idx_xstar : Indices of test items that are to be considered.
   *	\param[in] fbest :  An estimate of the best utility function value.
   *	\param[out] mei :  Maximum expected improvement for test items.
   */
    double maximum_expected_improvement(const VectorXd & theta_t, const VectorXd& theta_x, const double& sigma,
                                        const MatrixXd& t, const MatrixXd & x, const VectorXd& idx_global, const VectorXd& ind_t, const VectorXd& ind_x, MatrixXd tstar, int N, double fbest);

   /*! 
   * \brief Computes the expected value of information
   *
   *	Computes the expected value of information of including queries involving
   *	the pair test_pair. It  asssumes that:\n 
   *	-M is the number of training users + 1\n
   *	-t(M) = tstar\n
   *	-train_pairs(M) exists\n
   *
   *	\param[in] ind_t :  Indices of seen tasks.
   *	\param[in] ind_x :  Indices of seen items.
   *	\param[in] theta_x : The parameters for the items's covariance.
   *	\param[in] theta_t : The parameters for the users's covariance.
   *	\param[in] sigma : The noise parameter.
   *	\param[in] test_pair : The 2d vector representing the test query test_pair(1)>test_pair(2)
   *	\param[in] t :  The Matrix of User's features. 
   *	\param[in] x :  he Matrix of Item's features.
   *	\param[in] train_pairs : The Matrices of preferences pairs containing all the the preferences pairs of the users
   *	excepted the preferences pairs of the tested user.
   *	\param[in] idx_global : The unique global indices of the observed preferences.
   *	\param[in] tstar : Features of test user. 
   *	\param[in] idx_xstar : Indices of test items that are to be considered.
   *	\param[in] fbest :  An estimate of the best utility function value.
   *	\param[out] evoi: The expected value of information of asking the query
   *	test_pair(1)>test_pair(2).
   */
    double expected_voi(const VectorXd & theta_x, const VectorXd& theta_t, const double& sigma,
                        const MatrixXd& t, const MatrixXd & x, TypePair train_pairs, VectorXd& idx_global, VectorXd& ind_t, VectorXd& ind_x, MatrixXd test_pair, double fbest,double p_12);


   /*! 
   * \brief Computes The Elicitation Process for a given user
   *	
   *	Elicit preferences for a new user. It assumes that there is a
   *	TypePair Oracle passing in argument that knows
   *	the actual preference of the new user.
   *
   *	\param[in] theta_x : The parameters for the items's covariance.
   *	\param[in] theta_t : The parameters for the users's covariance.
   *	\param[in] sigma : The noise parameter.
   *	\param[in] idx_pairs : All the possible pairs.
   *	\param[in] train_pairs : The Matrices of preferences pairs containing all the the preferences pairs of the users
   *	excepted the preferences pairs of the tested user.
   *	\param[in] train_t :  The Matrix of training User's features. 
   *	\param[in] x :  he Matrix of Item's features.
   *	\param[in] train_pairs : The Matrices of preferences pairs containing all the the preferences pairs of the users
   *	excepted the preferences pairs of the tested user.
   *	\param[in] idx_global : The unique global indices of the observed preferences.
   *	\param[in] test_t : Features of test user. 
   *	\param[in] Maxiter : Maximum number of iterations (queries) to perform
   *	\param[in] Oracle :  The Matrices of preferences pairs containing all the the preferences pairs of the users
   *	including the preferences pairs of the tested user.
   *	\param[in] F : the targets Y: gaussian mean 0 cov K
   *	we sample this multivariate gaussian distribution N(mu,K)
   *	\param[out] loss : Vector of losses (one elemement per query)
   */
    void Elicit(const VectorXd & theta_x, const VectorXd& theta_t, const double& sigma, const MatrixXd& train_t, const MatrixXd &x, TypePair & train_pairs
                , const MatrixXd & test_t, int test_user_idx, MatrixXd  idx_pairs, int  Maxiter, const TypePair & Oracle, MatrixXd& F);



   /*! 
   * \brief Computes The Prediction Process for a given user
   *	
   *	Makes predictions on a new user by returning the test users' predicted utility at all items
   *
   *	\param[in] ind_t :  Indices of seen tasks.
   *	\param[in] ind_x :  Indices of seen items.
   *	\param[in] theta_x : The parameters for the items's covariance.
   *	\param[in] theta_t : The parameters for the users's covariance.
   *	\param[in] sigma : The noise parameter.
   *	\param[in] idx_pairs : All the possible pairs.
   *	\param[in] train_pairs : The Matrices of preferences pairs containing all the the preferences pairs of the users
   *	excepted the preferences pairs of the tested user.
   *	\param[in] train_t :  The Matrix of training User's features. 
   *	\param[in] x :  he Matrix of Item's features.
   *	\param[in] train_pairs : The Matrices of preferences pairs containing all the the preferences pairs of the users
   *	excepted the preferences pairs of the tested user.
   *	\param[in] idx_global : The unique global indices of the observed preferences.
   *	\param[in] test_t : Features of test user. 
   *	\param[in] idx_global: The unique global indices of the observed preferences
   *	\param[in] idx_global_1: The global indices of the first objects in the preferences
   *	\param[in] idx_global_2: idx_global_2: The gobal indices of the second objects in the preferences
   *	\param[in] ftrue: The true value of utility function of test user (used for evaluation)
   *	\param[in] ytrue: Binary vector indicating if the corresponding 
   *	item comparisons hold for the test user (used for evaluation)
   *	\param[out] fstar: the test users' predicted utility at all items
   */
    void Make_Predictions_New_User(const VectorXd & theta_x, const VectorXd& theta_t, double& sigma, const MatrixXd& train_t, const MatrixXd &x, const TypePair & train_pairs,
                                   const VectorXd & idx_global, const VectorXd& idx_global_1, const VectorXd& idx_global_2,
                                   const VectorXd& ind_t, const VectorXd& ind_x, const MatrixXd & test_t, const MatrixXd& idx_pairs, const VectorXd& ftrue, const VectorXd& ytrue);



   /*! 
   * \brief  Set the mean and variance of the predictive distribution of a gppe model
   *	
   *	Mustar and varstar are respectively set as the mean and the covariance of
   *	the predictive distribution.
   *
   *	\param[in] N : Number of Item's features.
   *	\param[in] idx_pairs : All the possible pairs.
   *	\param[in] t :  The Matrix of training User's features. 
   *	\param[in] idx_global : The unique global indices of the observed preferences.
   *	\param[in] test_t : Features of test user. 
   *	
   */
    void Predictive_Utility_Distribution(MatrixXd t, MatrixXd test_t, int N, VectorXd idx_global);


   /*! 
   * \brief Computes the probability p(x(test_pair(1),:) > x(test_pair(2),:) for a
   *	new user tstar using the parameters given by the Laplace method
   *
   *	Set p which is the probability that p(x(test_pair(1),:) > x(test_pair(2),:)\n
   *	and mustar which is the mean of the predictive distribution.	
   *
   *	\param[in] ind_t :  Indices of seen tasks.
   *	\param[in] ind_x :  Indices of seen items.
   *	\param[in] sigma : The noise parameter.
   *	\param[in] t :  The Matrix of training User's features. 
   *	\param[in] x :  he Matrix of Item's features.
   *	\param[in] idx_global : The unique global indices of the observed preferences.
   *	\param[in] tstar : Features of test user. 
   *	\param[in] test_pair = The test pair preference
   */
    void Predict_CGppe_Laplace(double sigma, MatrixXd t, MatrixXd x, VectorXd idx_global, VectorXd ind_t, VectorXd ind_x,
                              MatrixXd tstar, MatrixXd test_pair);


   /*! 
   * \brief Approximates the posterior distribution of the gppe model with the Laplace method
   *
   *	This procedure set the members of the CCPPE class which are : \n
   *	f: The mode of the posterior\n
   *	Kx: The covariance matrix on item space\n
   *	Kinv: The inverse covariance of the full system\n
   *	W: The matrix of negative second derivatives (wrt f) of the conditional likelihood\n
   *	L: chol(W + Kinv)' 
   *	llt : The Eigen Matrix Class which computes the Cholesky Factorisation
   *	
   *	\param[in] M :  The number of Users.
   *	\param[in] N :  The number of Items.
   *	\param[in] ind_t :  scalar index of the current user.
   *	\param[in] ind_x :  scalar index of the current item.
   *	\param[in] theta_x : The parameters for the items's covariance.
   *	\param[in] theta_t : The parameters for the users's covariance.
   *	\param[in] sigma : The noise parameter.
   *	\param[in] all_pairs : The Matrices of preferences pairs containing all the the preferences pairs of the users
   *	\param[in] t :  The Matrix of User's features. 
   *	\param[in] x :  he Matrix of Item's features.
   *	\param[in] idx_global : The unique global indices of the observed preferences.
   *	\param[in] idx_global_1: The global indices of the first objects in the preferences
   *	\param[in] idx_global_2: idx_global_2: The gobal indices of the second objects in the preferences
   */
    void Approx_CGppe_Laplace(const VectorXd& theta_x, const VectorXd& theta_t,
                             const double& sigma, const MatrixXd& t, const MatrixXd& x, const TypePair &all_pairs, const VectorXd &idx_global, const VectorXd& idx_global_1, const VectorXd& idx_global_2, const VectorXd& ind_t, const VectorXd& ind_x, int M, int N);


   /*! 
   * \brief Computes the conditional log-likelihood p(D | f,theta) of a gppe model
   *	
   *	\param[in] M :  The number of Users.
   *	\param[in] N :  The number of Items.
   *	\param[in] sigma : The noise parameter.
   *	\param[in] all_pairs : The Matrices of preferences pairs containing all the the preferences pairs of the users
   *	\param[in] idx_global_1: The global indices of the first objects in the preferences
   *	\param[in] idx_global_2: idx_global_2: The gobal indices of the second objects in the preferences

   *	\param[out] log_likelihood : the conditional log-likelihood p(D | f,theta)
   */
    double log_likelihood(double sigma, TypePair all_pairs, VectorXd idx_global_1, VectorXd idx_global_2, int M, int N);


   /*! 
   * \brief Computes the  derivatives of the conditional likelihood wrt f
   *	
   *	\param[in] M :  The number of Users.
   *	\param[in] N :  The number of Items.
   *	\param[in] sigma : The noise parameter.
   *	\param[in] all_pairs : The Matrices of preferences pairs containing all the the preferences pairs of the users
   *	\param[in] idx_global_1: The global indices of the first objects in the preferences
   *	\param[in] idx_global_2: idx_global_2: The gobal indices of the second objects in the preferences
   *	\param[out] deriv_loglike : the conditional log-likelihood p(D | f,theta)
   */
    VectorXd deriv_log_likelihood_CGppe_fast(double sigma, const TypePair & all_pairs, VectorXd idx_global_1, VectorXd idx_global_2, int M, int N);


   /*! 
   * \brief Computes the second derivatives of the conditional likelihood wrt f
   *	
   *	\param[in] M :  The number of Users.
   *	\param[in] N :  The number of Items.
   *	\param[in] sigma : The noise parameter.
   *	\param[in] all_pairs : The Matrices of preferences pairs containing all the the preferences pairs of the users
   *	\param[in] idx_global_1: The global indices of the first objects in the preferences
   *	\param[in] idx_global_2: idx_global_2: The gobal indices of the second objects in the preferences
   *	\param[out] Deriv2 : the second derivatives of the conditional likelihood wrt f
   */
    MatrixXd deriv2_log_likelihood_CGppe_fast(double sigma, TypePair all_pairs, VectorXd idx_global_1, VectorXd idx_global_2, int M, int N);


   /*! 
   * \brief Computes the expected best-so-far of the member f to get an expected loss
   *	
   *	\param[in] N :  The number of Items.
   *	\param[out] fbest : the expected best-so-far of the member f
   */
    double get_fbest(int N);
    
    
    //Members
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