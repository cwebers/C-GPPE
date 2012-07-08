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
#include "Tool.h"

class Covfunc
{
public :
    Covfunc();
    Covfunc(VectorXd theta);
    Covfunc(const Covfunc & c);
    ~Covfunc();
    VectorXd GetTheta();
    void SetTheta(VectorXd t);
    virtual int GetThetaLength() = 0;
    virtual double Evaluate(VectorXd t1, VectorXd t2) = 0;
    virtual double grad(VectorXd p, VectorXd q, int z)=0;
    MatrixXd ComputeGrandMatrix(MatrixXd fea);
    MatrixXd Compute(MatrixXd p, MatrixXd q);
    MatrixXd Computegrad(MatrixXd fea, int z);
    void add(int num);
protected :
    VectorXd Theta;
    //int a,b,c;
};

class CovSEard : public Covfunc
{
public :
    CovSEard(): Covfunc() {};
    CovSEard(VectorXd t): Covfunc(t) {};
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

class CovLINard : public Covfunc
{
public :
    CovLINard(): Covfunc() {};
    CovLINard(VectorXd t): Covfunc(t) {};
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

class CovNoise : public Covfunc
{
public :
    CovNoise(): Covfunc() {};
    CovNoise(VectorXd t): Covfunc(t) {};
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

class CovSEiso : public Covfunc
{
public :
    CovSEiso(): Covfunc() {};
    CovSEiso(VectorXd t): Covfunc(t) {};
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

class CovSum : public Covfunc
{
public :
    CovSum(Covfunc *ptr1, Covfunc *ptr2): Covfunc()
    {
        p_a = ptr1;
        p_b = ptr2;
    };
    CovSum(Covfunc *ptr1, Covfunc *ptr2, VectorXd t): Covfunc(t)
    {
        p_a = ptr1;
        p_b = ptr2;
    };
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
