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
#ifndef __Covfunc_H__
#define __Covfunc_H__
#include <Eigen/Dense>
#include "dlib/optimization.h"
#include <iostream>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
using namespace std;

/*! \namespace Eigen
 * \brief Namespace gathering the containers of the Eigen Library and its Optimizers.
 *
 *	Eigen is a C++ template library for linear algebra: matrices, vectors, numerical solvers, and related algorithms.
 *	For more information, please vite the Eigen documentation website : http://eigen.tuxfamily.org/dox/
 *	
 */ 
 using namespace Eigen;

 /*! \namespace Dlib
 * \brief Namespace gathering the containers of the Dlib Library and its Optimizers.
 *
 *	Dlib is a general purpose cross-platform C++ library designed using contract programming and modern C++ techniques.
 *	For more information, please visit the Dlib documentation websites : http://dlib.net/
 */
using namespace dlib;


#include <iostream>
#include <time.h>
#include <sstream>



/*! \class VectorXd
   * \brief Dynamic Vector of double
   *
   *	This Class implemented by Eigen provide a Dynamic Vector of Double
   *	possessing a large variety of methods	
   */
   
   /*! \class MatrixXd
   * \brief Dynamic Matrix of double
   *
   *	This Class implemented by Eigen provide a Dynamic Matrix of Double
   *	possessing a large variety of methods	
   */


/*! \class Covfunc
   * \brief Virtual class for the Covariance
   *
   *  This Class gather the major methods that all the Covariances classes share together
   * and the Parameter Theta as a member of the Class.
   */
class Covfunc
{
public :

/*!
     *  \brief Constructor by default
     *
     *  Constructor of the Covfunc Class
     *
     */
    Covfunc();
    
    /*!
     *  \brief Instanciate with the parameter Theta.
     *
     *  Constructor of the Covfunc Class.
     *	\param VectorXd Theta
     */
    Covfunc(const VectorXd& theta);
    
        /*!
     *  \brief Constructor by recopy
     *
     *  Instanciate an object of the class Covfunc by copying an
     *	other object of the class.
     *	\param Covfunc c
     */
    Covfunc(const Covfunc & c);
    
        /*!
     *  \brief Destructor of the Covfunc Class
     */
    ~Covfunc();
    
        /*!
     *  \brief Accessor for Theta.
     *
     *  Return the member Theta.
     */
    VectorXd GetTheta();
    
     /*!
     *  \brief Mutator for Theta.
     *
     *  Set the member Theta with the given parameter t.
     * \param VectorXd t
     */
    void SetTheta(const VectorXd t);
   
     /*!
     *  \brief Return the length of Theta following a code
     *
     *  Return an integer code following the type of the Covariance.
     *	If the parameter code is negative, then the length of theta is
     *	abs(code)+lenght of the vector we are dealing with.
     *	If it's zero, the length is which of the vector we are dealing with.
     *  If it's superior to zero, the code is the actual length of Theta.
     */
    virtual int GetThetaLength() = 0;
    
     /*!
     *  \brief Evaluate the Covariance between 2 vectors.
     *
     *  Calculate the Covariance between 2 vectors following theta and the
     *	type of Covariance class used.
     *	\param VectorXd t1, VectorXd t2
     */
    virtual double Evaluate(VectorXd t1, VectorXd t2) = 0;
    
         /*!
     *  \brief Evaluate the gradient of the Covariance between 2 vectors.
     *
     *  Calculate the gradient of the Covariance between 2 vectors, following theta, 
     *  the type of Covariance class used and the index of the parameter in theta
     * which you want to derivate with.
     *	\param VectorXd p, VectorXd q, int z
     */
    virtual double grad(VectorXd p, VectorXd q, int z)=0;
      
         /*!
     *  \brief Return the Matrix of covariances
     *	
     * Return a Matrix Which contains the covariances between all
     * the Vectors inside the matrix fea.
     *	\param MatrixXd fea
     */
    MatrixXd ComputeGrandMatrix( MatrixXd fea);
  
            /*!
     *  \brief Return the Matrix of covariances
     *	
     * Return a Matrix Which contains the covariances between all
     * the Vectors inside the matrix p and the vectors inside q
     *	\param MatrixXd p, MatrixXd q
     */
    MatrixXd Compute(MatrixXd p, MatrixXd q);
      
         /*!
     *  \brief Return the gradient Matrix of covariances
     *	
     * Return a Matrix Which contains the gradient of the covariances between all
     * the Vectors inside the matrix fea. Derivate following the index z
     *	\param MatrixXd fea, int z
     */
    MatrixXd Computegrad(MatrixXd fea, int z);
protected :
    VectorXd Theta;
};


/*! \class CovSEard
   * \brief Class which computes the CovSEard
   *
   * This class computes the Squared Exponential covariance 
   * function with Automatic Relevance Determination
   */
class CovSEard : public Covfunc
{
public :
         /*!
     *  \brief Constructor by default
     */
    CovSEard(): Covfunc() {};
    
    /*!
     *  \brief Instanciate with the parameter Theta.
     *
     *  Constructor of the CovSEard Class which initialize Theta with the given parameter.
     *	\param VectorXd Theta
     */
    CovSEard(VectorXd t): Covfunc(t) {};
    
            /*!
     *  \brief Destructor of the CovSEard Class
     */
    ~CovSEard() {};

    virtual double Evaluate(VectorXd p, VectorXd q)
    {
        double res = 0;
        double sum = 0;
        double sf2 = exp(2 * Theta(Theta.rows() - 1));
        for (int i = 0;i < p.rows();i++)
        {
            sum += pow((p(i) - q(i)), 2) / exp(2 * Theta(i));
        }
        res = sf2 * exp(-0.5 * sum);


        return res;
    }
    
    virtual double grad(VectorXd p, VectorXd q, int z)
    {
        double res = 0;
        double deriv=pow(exp(Theta(z)),2);
        double prodsup;
        double sum = 0;
        double sf2 = exp(2 * Theta(Theta.rows() - 1));
        if(z==Theta.rows()-1)
        	return 2*Evaluate(p,q)/sf2;
        	
        for (int i = 0;i < p.rows();i++)
        {
            sum += pow((p(i) - q(i)), 2) / exp(2 * Theta(i));
        }
        prodsup=pow((p(z) - q(z)), 2)/pow(deriv,2);
        res = prodsup*sf2 * exp(-0.5 * sum);
        return res;
    }
        
    virtual int GetThetaLength()
    {
        return -1;
    }
};

/*! \class CovLINard
   * \brief Class which computes the CovLINard
   *
   *  This class implements the Linear covariance function with Automatic Relevance Determination (ARD).
   */
class CovLINard : public Covfunc
{
public :

         /*!
     *  \brief Constructor by default
     */
    CovLINard(): Covfunc() {};
    
    /*!
     *  \brief Instanciate with the parameter Theta.
     *
     *  Constructor of the CovLINars Class which initialize Theta with the given parameter.
     *	\param VectorXd Theta
     */
    CovLINard(VectorXd t): Covfunc(t) {};
    
            /*!
     *  \brief Destructor of the CovLINard Class
     */
    ~CovLINard() {};
    virtual double Evaluate(VectorXd p, VectorXd q)
    {
        double res = 0;
        for (int i = 0;i < p.rows();i++)
        {
            res += (p(i) * q(i)) / (exp(2 * Theta(i)));
        }

        return res;
    }
    
    virtual double grad(VectorXd p, VectorXd q, int z)
    {
    	return -1;
    }
    
    virtual int GetThetaLength()
    {
        return 0;
    }
};

/*! \class CovNoise
   * \brief Class which computes the CovNoise
   *
   *  This class implements the Covariance Noise.
   */
class CovNoise : public Covfunc
{
public :

         /*!
     *  \brief Constructor by default
     */
    CovNoise(): Covfunc() {};
    
        /*!
     *  \brief Instanciate with the parameter Theta.
     *
     *  Constructor of the CovNoise Class which initialize Theta with a given parameter.
     *	\param VectorXd Theta
     */
    CovNoise(VectorXd t): Covfunc(t) {};
    
            /*!
     *  \brief Destructor of the CovNoise Class
     */
    ~CovNoise() {};
    virtual double Evaluate(VectorXd p, VectorXd q)
    {
        double res = 0;
        if (p == q)
            return exp(2*Theta(0));
        else
            return 0;
    }
    
    virtual double grad(VectorXd p, VectorXd q, int z)
    {
    	return -1;
    }
    
    virtual int GetThetaLength()
    {
        return 1;
    }
};

/*! \class CovSEiso
   * \brief Class which computes the CovSEiso
   *
   *  This class implements the Squared Exponential covariance function with isotropic distance measure.
   */
class CovSEiso : public Covfunc
{
public :

         /*!
     *  \brief Constructor by default
     */
    CovSEiso(): Covfunc() {};
    
        /*!
     *  \brief Instanciate with the parameter Theta.
     *
     *  Constructor of the CovSEiso Class which initialize Theta with a given parameter.
     *	\param VectorXd Theta
     */
    CovSEiso(VectorXd t): Covfunc(t) {};
    
            /*!
     *  \brief Destructor of the CovSEiso Class
     */
    ~CovSEiso() {};

    virtual double Evaluate(VectorXd p, VectorXd q)
    {
        double res = 0;
        double sum = 0;
        double sf2 = exp(2 * Theta(1));
        for (int i = 0;i < p.rows();i++)
        {
            sum += pow((p(i) - q(i)), 2);
        }
        res = sf2 * exp(-0.5 * sum / exp(2 * Theta(0)));


        return res;
    }
    
    
    virtual double grad(VectorXd p, VectorXd q, int z)
    {
    	return -1;
    }
    
    virtual int GetThetaLength()
    {
        return 2;
    }
};

/*! \class CovSum
   * \brief Class which computes the CovSum
   *
   *  This class takes two types of covariances and Evaluates as the sum of theses covariances Evaluates methods
   */
class CovSum : public Covfunc
{
public :


         /*!
     *  \brief Constructor by default
     *
     *Take two 2 Covariances type as members
     */
    CovSum(Covfunc *ptr1, Covfunc *ptr2): Covfunc()
    {
        p_a = ptr1;
        p_b = ptr2;
    };
    
    /*!
     *  \brief Instanciate with the parameter Theta.
     *
     *  Constructor of the CovSum Class Which initialize The parameters of the 2 covariances
     *	members. The parameters given contains the two Thetas concatained.
     *	\param VectorXd Theta
     */
    CovSum(Covfunc *ptr1, Covfunc *ptr2, VectorXd t): Covfunc(t)
    {
        p_a = ptr1;
        p_b = ptr2;
    };
    
                /*!
     *  \brief Destructor of the CovSum Class
     */
    ~CovSum() {};

    virtual double Evaluate(VectorXd p, VectorXd q)
    {
        switch (p_a->GetThetaLength())
        {
        case -1 :
            p_a->SetTheta(Theta.block(0, 0, p.rows() + 1, 1));
            p_b->SetTheta(Theta.block(p.rows() + 1, 0, GetThetaLength() - p.rows() - 1, 1));
            break;

        case 0 :
            p_a->SetTheta(Theta.block(0, 0, p.rows(), 1));
            p_b->SetTheta(Theta.block(p.rows(), 0, GetThetaLength() - p.rows(), 1));
            break;

        default:
            p_a->SetTheta(Theta.block(0, 0, p_a->GetThetaLength(), 1));
            p_b->SetTheta(Theta.block(p_a->GetThetaLength(), 0, GetThetaLength() - p_a->GetThetaLength(), 1));
        }

        return p_a->Evaluate(p, q) + p_b->Evaluate(p, q);
    }
    
    virtual double grad(VectorXd p, VectorXd q, int z)
    {
    	return -1;
    }

    virtual int GetThetaLength()
    {
        return Theta.rows();
    }
protected :
    Covfunc *p_a;
    Covfunc *p_b;

};
#endif
