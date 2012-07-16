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
#include <iostream>
#include <string>
#include "CGppe.h"


/*! \class CLearner
   * \brief Class for the hyperparameter's Optimization
   *
   *	This class computes the negative_marginal_log_likelihood and its gradient.
   *	It is used in an Optimizer
   */
class CLearner
{
public :

/*! 
   * \brief Constructor by default
   *
   */
    CLearner();
    /*!    
     * \brief Constructor with instanciation of the members
     *
     *	\param[in] ind_t :  scalar index of the current user
     *	\param[in] ind_x :  scalar index of the current item
     *	\param[in] all_pairs : TypePair of M elements. Each element is a O_m x 2 matrix 
     *	where O_m is the number of preferences observed for the corresponding
     *	user. Each row all_pairs(m) contains a preference relation 
     *	of the form all_pairs(m)(0) > all_pairs(m)(1) 
     *	\param[in] M : The number of users
     *	\param[in] N : The number of items 
     *	\param[in] idx_global : The unique global indices of the observed preferences
     *	\param[in] idx_global_1 : The global indices of the first objects in the preferences
     *	\param[in] idx_global_2 : The gobal indices of the second objects in the preferences
   */
    CLearner(Covfunc *Covt, Covfunc *Covx,
          MatrixXd T, MatrixXd X, TypePair All_pairs, VectorXd Idx_global, VectorXd Idx_global_1,
          VectorXd Idx_global_2, VectorXd Ind_t, VectorXd Ind_x, int  m, int n);// Default Contructor

/*! 
   * \brief Destructor
   *
   */
    ~CLearner();


    /*!    
     * \brief Returns the negative marginal log likelihood
     *
     *	Returns the negative marginal log likelihood given this instanciates parameters and
     *	the input theta.
     *	\param[in] dltheta : The Hyperparameters theta_t, theta_x and logsigma concatenated into a single Vector
     *	\param[out] : The negative marginal log likelihood.  
   */
    const double negative_marginal_log_likelihood(const column_vector & dltheta);
    
        /*!    
     * \brief Returns the gradient of the negative marginal log likelihood
     *
     *	Returns the gradient of the negative marginal log likelihood given this instanciates parameters and
     *	the input theta.
     *	\param[in] dltheta : The Hyperparameters theta_t, theta_x and logsigma concatenated into a single Vector
     *	\param[out] : The gradient of the hyperparameters.  
   */
    column_vector gradient_negative_marginal_loglikelihood(const column_vector & point);
    
            /*!    
     * \brief Accesor for the min of the negative marginal log likelihood 
     *
     *	Returns the min of the negative marginal log likelihood once the optimization has been done.
     *	\param[out] : The min of the negative marginal log likelihood .  
   */
	double Getmin();

    
  
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
    double min;


};

/*! \class CNLL_Function
   * \brief Wrapper for negative log likelihood
   *
   *	This object is an example of what is known as a "function object" in C++.
   *	It is simply an object with an overloaded operator(). This means it can
   *	be used in a way that is similar to a normal C function. The interesting
   *	thing about this sort of function is that it can have state.Class used to overload
   *	the negative_marginal_loglikelihood of the  CLearner CLass.
   */
class CNLL_Function
{
public:
    /*!    
   * \brief Constructor with instanciation of a CLearner Object
   */       
    CNLL_Function ( const CLearner & learner) : _learner(learner)
    {
    }

        /*!    
   * \brief Overload the negative_marginal_loglikelihood's function of the
   *	CLearner member.
   */
    double operator() ( const column_vector & arg) const
    {
        return ((CLearner)_learner).negative_marginal_log_likelihood(arg);
    }
    
private:
    CLearner _learner;
};




/*! \class CGradNLL_Function
   * \brief Wrapper for grad negative log likelihood
   *
   *	This object is an example of what is known as a "function object" in C++.
   *	It is simply an object with an overloaded operator(). This means it can
   *	be used in a way that is similar to a normal C function. The interesting
   *	thing about this sort of function is that it can have state. Class used to overload
   *	the gradient_negative_marginal_loglikelihood of the  CLearner CLass.
   */
class CGradNLL_Function
{
public:
    /*!    
   * \brief Constructor with instanciation of a CLearner Object
   */   
    CGradNLL_Function ( const CLearner & learner) :_learner(learner)
    {
    }
        /*!    
   * \brief Overload the gradient_negative_marginal_loglikelihood's function of the
   *	CLearner member.
   */   
    column_vector operator() ( const column_vector & arg) const
    {
       return ((CLearner)_learner).gradient_negative_marginal_loglikelihood(arg);
    }
    
private:
    CLearner _learner;
};
#endif