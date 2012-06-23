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
#include "Covfunc.h"

Covfunc::Covfunc() 
{
	Theta= VectorXd(3);
	//a=5;
	//b=2;
	//c=0;
}

Covfunc::Covfunc(const Covfunc & c)
{
	VectorXd Theta;
	Theta=Theta=c.Theta;
}

Covfunc::Covfunc( VectorXd theta) 
{
	Theta = theta;
}

Covfunc::~Covfunc() 
{

}

VectorXd Covfunc::GetTheta() 
{ 
	return Theta; 
}


void Covfunc::SetTheta(VectorXd t) 
{
	Theta=t;
}

/*void Covfunc::add(int num) 
{
	c=a+b+num;
	cout<<c<<endl;
}*/


MatrixXd Covfunc::ComputeGrandMatrix(MatrixXd fea)
{
	MatrixXd cov(fea.rows(),fea.rows());
	for(int i=0;i<fea.rows();i++)
	{
		for(int j=0;j<fea.rows();j++)
		{
			cov(i,j)=Evaluate(fea.row(i),fea.row(j));
		}
	}
	return cov;
}

MatrixXd Covfunc::Compute(MatrixXd p, MatrixXd q)
{
	MatrixXd cov(p.rows(),p.rows());
	for(int i=0;i<p.rows();i++)
	{
		for(int j=0;j<p.rows();j++)
		{
			cov(i,j)=Evaluate(p.row(i),q.row(j));
		}
	}
	return cov;
}

double normcdf(double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return 0.5*(1.0 + sign*y);
}

VectorXd normcdf(VectorXd x)
{
	for(int i=0;i<x.rows();i++)
	{
		x(i)=normcdf(x(i));
	}
	return x;
}


