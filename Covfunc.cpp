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
	MatrixXd cov=MatrixXd::Zero(fea.rows(),fea.rows());
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
	MatrixXd cov=MatrixXd::Zero(p.rows(),p.rows());
	for(int i=0;i<p.rows();i++)
	{
		for(int j=0;j<p.rows();j++)
		{
			cov(i,j)=Evaluate(p.row(i),q.row(j));
		}
	}
	return cov;
}

