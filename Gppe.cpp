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
#include "Gppe.h"

Gppe::Gppe()
{

}



Gppe::Gppe(const Gppe & g)
{
	f=g.f;
	Kx=g.Kx;
	Kinv=g.Kinv;
	W=g.W;
	L=g.L;
}

Gppe::Gppe(VectorXd fnew,MatrixXd Kxnew,MatrixXd Kinvnew,MatrixXd Wnew,MatrixXd Lnew) 
{
	f=fnew;
	Kx=Kxnew;
	Kinv=Kinvnew;
	W=Wnew;
	L=Lnew;

}

Gppe::~Gppe() 
{

}







void Gppe::Approx_Gppe_Laplace(Covfunc *Covfunc_t,Covfunc *Covfunc_x,
			VectorXd theta_x,VectorXd theta_t, double sigma,MatrixXd t,MatrixXd x,TypePair all_pairs,
			VectorXd idx_global,VectorXd idx_global_1,VectorXd idx_global_2, 
			VectorXd ind_t,VectorXd ind_x,int M,int N)
{
	//Parameters function initialization
	double eps=10E-6, psi_new,psi_old;
	int n=M*N;
	sigma=exp(sigma);
	VectorXd f(n);
	VectorXd fvis(idx_global.rows());

	Covfunc_t->SetTheta(theta_t);
	Covfunc_x->SetTheta(theta_x);

	MatrixXd Kt=Covfunc_t->ComputeGrandMatrix(t);
	Kx=Covfunc_x->ComputeGrandMatrix(x);

	MatrixXd K= ComputeMat(Kt,ind_t,ind_t).array()*ComputeMat(Kx,ind_x,ind_x).array();
	
	double loglike =log_likelihood(f, sigma, all_pairs, idx_global_1, idx_global_2, N);
	Kinv=K.inverse();
	psi_new = loglike - 0.5 * fvis.transpose()*Kinv*fvis; 
	psi_old = -1000;
	while((psi_new-psi_old)>eps)
	{
		psi_old=psi_new;
		
	}

}

MatrixXd Gppe::ComputeMat(MatrixXd mat,VectorXd t1, VectorXd t2)
{
	MatrixXd res(t1.rows(),t2.rows());
	for(int i=0;i<t1.rows();i++)
	{
		for(int j=0;j<t2.rows();j++)
		{
			res(i,j)=mat(t1(i),t2(j));
		}
	}
	return res;
}

VectorXd Gppe::ComputeVec(VectorXd vec,VectorXd t1)
{
	VectorXd res(t1.rows());
	for(int i=0;i<t1.rows();i++)
	{
		res(i)=vec(t1(i));
	}
	return res;
}

VectorXd Gppe::ind2global(VectorXd vec,int j,int N)
{
	VectorXd res;
	res= (j-1)*N +vec.array();
	
	return res;
}


double Gppe::log_likelihood(VectorXd f,double sigma, TypePair all_pairs,VectorXd idx_global_1,VectorXd idx_global_2, int N)
{
	MatrixXd pairs;
	VectorXd idx_1,idx_2,z;
	double loglike=0;
	int M=all_pairs.rows();
	for(int j=0;j<M;j++)
		{	
			if(all_pairs(j).rows()==0)
				continue;
			pairs=all_pairs(j);
			idx_1=ind2global(pairs.col(1),j,N);
			idx_2=ind2global(pairs.col(2),j,N);
			z=(ComputeVec(f,idx_1)-ComputeVec(f,idx_2))/sigma;
			z=normcdf(z);
			loglike+= log(z.array()).sum();	
		}
	return loglike;
}
