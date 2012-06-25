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


VectorXd Gppe::Getf()
{
	return f;
}

MatrixXd Gppe::GetW()
{
	return W;
}

MatrixXd Gppe::GetL()
{
	return L;
}

MatrixXd Gppe::GetKinv()
{
	return Kinv;
}

MatrixXd Gppe::GetKx()
{
	return Kx;
}
void Gppe::Approx_Gppe_Laplace(Covfunc *Covfunc_t,Covfunc *Covfunc_x,
			VectorXd theta_x,VectorXd theta_t, double sigma,MatrixXd t,MatrixXd x,TypePair all_pairs,
			VectorXd idx_global,VectorXd idx_global_1,VectorXd idx_global_2, 
			VectorXd ind_t,VectorXd ind_x,int M,int N)
{
	//Parameters function initialization
	double eps=10E-6, psi_new,psi_old;
	int n=M*N;
	//sigma=exp(sigma);
	VectorXd f(n);
	VectorXd fvis(idx_global.rows());
	VectorXd deriv;
	LLT<MatrixXd> llt;

	Covfunc_t->SetTheta(theta_t);
	Covfunc_x->SetTheta(theta_x);

	MatrixXd Kt=Covfunc_t->ComputeGrandMatrix(t);
	Kx=Covfunc_x->ComputeGrandMatrix(x);

	MatrixXd K= GetMat(Kt,ind_t,ind_t).array()*GetMat(Kx,ind_x,ind_x).array();
	double loglike =log_likelihood(f, sigma, all_pairs, idx_global_1, idx_global_2, N);
	Kinv=K.inverse();
	psi_new = loglike - 0.5 * fvis.transpose()*Kinv*fvis; 
	psi_old = -10E6;
	while((psi_new-psi_old)>eps)
	{
		psi_old=psi_new;
		deriv=deriv_log_likelihood_gppe_fast(f, sigma, all_pairs, idx_global_1, idx_global_2, M, N);
		W=-deriv2_log_likelihood_gppe_fast(f, sigma, all_pairs, idx_global_1, idx_global_2, M, N);
		cout<<W<<endl;
		W=GetMat(W,idx_global,idx_global);
		llt.compute(W+Kinv);
		L = llt.matrixL(); //no need to extract the triangular matrix here
		fvis=llt.solve(GetVec(deriv,idx_global)+W*fvis);
		
		for(int w=0;w<idx_global.rows();w++)
		{
			f(idx_global(w))=fvis(w);	
		}
		loglike =log_likelihood(f, sigma, all_pairs, idx_global_1, idx_global_2, N);
		psi_new = loglike - 0.5 * fvis.transpose()*Kinv*fvis; 
	}

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
			idx_1=ind2global(pairs.col(0),j,N);
			idx_2=ind2global(pairs.col(1),j,N);
			z=(GetVec(f,idx_1)-GetVec(f,idx_2))/sigma;
			z=normcdf(z);
			loglike+= log(z.array()).sum();	
		}
	return loglike;
}

VectorXd Gppe::deriv_log_likelihood_gppe_fast(VectorXd f,double sigma,const  TypePair& all_pairs, VectorXd idx_global_1, VectorXd idx_global_2, int M, int N)
{
	VectorXd deriv_loglike, z, cdf_val, pdf_val,val;
	// test if idx vectors are empty ?
	M=all_pairs.rows();
	int n=M*N;
	z=(GetVec(f,idx_global_1)-GetVec(f,idx_global_2))/sigma;
	cdf_val= normcdf(z);
	pdf_val= normpdf(z);
	val= (1./sigma) * (pdf_val.array()/cdf_val.array());
	return Get_Cumulative_Val(idx_global_1, val, n)- Get_Cumulative_Val(idx_global_2, val, n);
}

MatrixXd Gppe::deriv2_log_likelihood_gppe_fast(VectorXd f,double sigma, TypePair all_pairs, VectorXd idx_global_1, VectorXd idx_global_2, int M, int N)
{
	VectorXd deriv_loglike, z, cdf_val, pdf_val,val,ratio,all_diag_idx,ind,ind_trans;
	M=all_pairs.rows();
	int n=M*N;
	VectorXd consec(n);
	MatrixXd Deriv2(n,n);
	for(int i=0;i<n;i++)
	{
		consec(i)=i;
	}
	all_diag_idx=sub2ind(n,n,consec,consec);
	z=(GetVec(f,idx_global_1)-GetVec(f,idx_global_2))/sigma;
	cdf_val= normcdf(z);
	pdf_val= normpdf(z);
	ratio=pdf_val.array()/cdf_val.array();
	val= -(1./pow(sigma,2))*(ratio.array()*(z+ratio).array());
	
	ind=sub2ind(n,n,idx_global_1,idx_global_2);
	GetMatGenIdx(Deriv2,ind,-val);
	ind_trans=sub2ind(n,n,idx_global_2,idx_global_1);
	GetMatGenIdx(Deriv2,ind_trans,-val);
	
	GetMatGenIdx(Deriv2,all_diag_idx,GetMatGenIdx(Deriv2,all_diag_idx)+Get_Cumulative_Val(idx_global_1, val, n));
	GetMatGenIdx(Deriv2,all_diag_idx,GetMatGenIdx(Deriv2,all_diag_idx)+Get_Cumulative_Val(idx_global_2, val, n));
	return Deriv2;
}
