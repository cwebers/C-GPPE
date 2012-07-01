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

Gppe::Gppe(Covfunc *Covfunc_t,Covfunc *Covfunc_x)
{
	covfunc_t=Covfunc_t;
	covfunc_x=Covfunc_x;
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

VectorXd Gppe::Getmustar()
{
	return mustar;
}

VectorXd Gppe::Getvarstar()
{
	return varstar;
}


double Gppe::get_fbest(int N)
{
	VectorXd ftest;
	double fbest;
	
	ftest=f.segment(f.rows()-N,N-1);
	fbest=ftest.maxCoeff();
	if(fbest!=fbest)
		fbest=f.maxCoeff();
	return fbest;
}


void Gppe::Elicit( const VectorXd & theta_x,const VectorXd& theta_t, const double& sigma,const MatrixXd& train_t,const MatrixXd &x,const TypePair & train_pairs
    , const MatrixXd & test_t, int test_user_idx, MatrixXd  idx_pairs,int  Maxiter)// ptr_query_func, ptr_loss_func)
{
	int N=x.rows();
	int Mtrain=train_t.rows();
	int M=Mtrain +1;
	int Npairs=idx_pairs.rows();
	//VectorXd is_selected(Npairs);
	Matrix<bool,Dynamic,1> is_selected(Npairs);
	is_selected.fill(false);
	VectorXd loss=VectorXd::Zero(Maxiter+1);
	VectorXd evoi(Npairs), ind_t, ind_x;
	VectorXd idx_global_1, idx_global_2, idx_global;
	compute_global_index(idx_global_1,idx_global_2,train_pairs, N);
	unique(idx_global,idx_global_1,idx_global_2);
	bool stop=false;
	double foo, val, idx_good;
	int count=0;
	MatrixXd t;
	
	
	t.resize(M,train_t.cols());
	t<<train_t,test_t;
	TypePair inter(M);
	inter<<train_pairs;//need to check this entry later	
	
	for(int iter=0;iter<=Maxiter;iter++)
	{	
	
		Approx_Gppe_Laplace( theta_x, theta_t, sigma,
    	t, x, train_pairs, idx_global, idx_global_1, idx_global_2, ind_t, ind_x, Mtrain, N);
	
	
		Predictive_Utility_Distribution(t, test_t,N,idx_global );
	
		std::ptrdiff_t idxbest;
		foo = mustar.maxCoeff(&idxbest);
		double fbest=get_fbest(N);
	
		for(int i=0;i<Npairs;i++)
		{
			if(is_selected(i))
			{
				evoi(i)=INT_MIN;
				continue;
			}
			VectorXd test_pair=idx_pairs.row(i);
			//evoi(i)=expected_voi();
		}
		
		std::ptrdiff_t query_idx;
		val = evoi.maxCoeff(&query_idx);
		idx_good=find(evoi,val);
		//evoi is a vector, so there isn't any multiple index argument
		//and no need for the Lgood check stuff
		is_selected(query_idx)=false;
		
		
		
		
		
		
		count++;
	}
	
	
	
	
	
}


void Gppe::Make_Predictions_New_User(const VectorXd & theta_x,const VectorXd& theta_t, const double& sigma,const MatrixXd& train_t,const MatrixXd &x,const TypePair & train_pairs,
			const VectorXd & idx_global,const VectorXd& idx_global_1,const VectorXd& idx_global_2, 
			const VectorXd& ind_t,const VectorXd& ind_x,const MatrixXd & test_t,const MatrixXd& idx_pairs,const VectorXd& ftrue, const VectorXd& ytrue)
{
	int N=x.rows();
	int Mtrain=train_t.rows();
    int Npairs= idx_pairs.rows();
    VectorXd fstar;
    //MatrixXd Fstar=MatrixXd::Zero(N,Npairs);
    //Fstar=SetNaN(Fstar);
    VectorXd pair;
    VectorXd P=VectorXd::Zero(Npairs);
    VectorXd ypred=VectorXd::Zero(Npairs);
    VectorXd sum=VectorXd::Zero(N);
    VectorXd count=VectorXd::Zero(N);
    Approx_Gppe_Laplace( theta_x, theta_t, sigma,
    train_t, x, train_pairs, idx_global, idx_global_1, idx_global_2, ind_t, ind_x, Mtrain, N);
    
    for(int i=0;i<Npairs;i++)
    {
    	pair=idx_pairs.row(i);
    	Predict_Gppe_Laplace(sigma, train_t,x,idx_global,ind_t, ind_x,
		 test_t, pair);
		 P(i)=p;
		 //Fstar(pair(0),i)=mustar(0);
		 //Fstar(pair(1),i)=mustar(1);
		 sum(pair(0))+=mustar(0);
		 count(pair(0))+=1;
		 sum(pair(1))+=mustar(1);
		 count(pair(1))+=1;
		 
    }  
    
    for(int i=0;i<P.rows();i++)
    {
    	if(P(i)>0.5)
    		ypred(i)=1;
    } 

//fstar=MyNaNMean(Fstar);
fstar=sum.array()/count.array();
    
cout<<endl<<endl<<"error =  "<<(GetDiff(ytrue,ypred)).sum()/ytrue.rows();
  // need for a plot function here ?  
}





void Gppe::Predictive_Utility_Distribution(MatrixXd t,MatrixXd tstar, int N, VectorXd idx_global)
{
	int Kt_ss=1;
	VectorXd idx_xstar(N);
	MatrixXd Kstar,Kx_star_star,Kx_star,Kss,Css,Kt_star;
	for(int i=0;i<N;i++)
	{
		idx_xstar(i)=i;
	}
	
	Kt_star=covfunc_t->Compute(t,tstar);
	Kx_star=GetMatRow(Kx,idx_xstar);//need to check for tranpose later?
	Kx_star_star=GetMat(Kx,idx_xstar, idx_xstar);
	
	Kstar=Kron(Kt_star,Kx_star);

	Kstar=GetMatRow(Kstar,idx_global);

	Kss=Kt_ss*Kx_star_star;
	mustar=Kstar.transpose()*Kinv*GetVec(f,idx_global);
	Css= Kss- Kstar.transpose()*W*llt.solve(Kinv*Kstar);
	varstar=Css.diagonal();
}


void Gppe::Predict_Gppe_Laplace(double sigma, MatrixXd t, MatrixXd x, VectorXd idx_global, VectorXd ind_t, VectorXd ind_x,
MatrixXd tstar, MatrixXd test_pair)
{
	int Kt_ss=1;
	double sigma_star, val;
	MatrixXd Kx_star,Kx_star_star,kstar,Kss,Css;
	MatrixXd Kt_star=covfunc_t->Compute(t,tstar);


	
	Kx_star = GetMatRow(Kx,test_pair.transpose()); //maybe need some transpose?
	cout<<Kt_star.cols()<<endl;
	
	Kx_star_star = GetMat(Kx,test_pair.transpose(),test_pair.transpose()); // test to test
	kstar = Kron(Kt_star, Kx_star);
	
	kstar =GetMatRow(kstar,idx_global);
	Kss = Kt_ss * Kx_star_star;


	mustar = kstar.transpose()*Kinv*GetVec(f,idx_global);
	Css    = Kss - kstar.transpose()*W*llt.solve(Kinv*kstar); 

	sigma_star = sqrt(Css(0,0) + Css(1,1) - 2*Css(0,1) + pow(sigma,2));
	cout<<sigma_star<<endl;
	val = ( mustar(0) - mustar(1) )/sigma_star;
	p   = normcdf(val);
}




void Gppe::Approx_Gppe_Laplace(const VectorXd & theta_x,const VectorXd& theta_t, const double& sigma,const MatrixXd& t,const MatrixXd &x,const TypePair & all_pairs,
			const VectorXd & idx_global,const VectorXd& idx_global_1,const VectorXd& idx_global_2, 
			const VectorXd& ind_t,const VectorXd& ind_x,int M,int N)
{
	//Parameters function initialization
	double eps=1E-6, psi_new,psi_old;
	M=all_pairs.rows();
	int n=M*N;
	 f=VectorXd::Zero(n);
	VectorXd fvis=VectorXd::Zero(idx_global.rows());
	VectorXd deriv;
	double loglike=0;

	covfunc_t->SetTheta(theta_t);
	covfunc_x->SetTheta(theta_x);

	MatrixXd Kt=covfunc_t->ComputeGrandMatrix(t);
	Kx=covfunc_x->ComputeGrandMatrix(x);

	MatrixXd K= GetMat(Kt,ind_t,ind_t).array()*GetMat(Kx,ind_x,ind_x).array();

	loglike = log_likelihood(f, sigma, all_pairs, idx_global_1, idx_global_2,M, N);
	Kinv=K.inverse();
	psi_new = loglike - 0.5 * fvis.transpose()*Kinv*fvis; 
	psi_old = INT_MIN;
	while((psi_new-psi_old)>eps)
	{
		psi_old=psi_new;
		deriv=deriv_log_likelihood_gppe_fast(f, sigma, all_pairs, idx_global_1, idx_global_2, M, N);
		W=-deriv2_log_likelihood_gppe_fast(f, sigma, all_pairs, idx_global_1, idx_global_2, M, N);
		W=GetMat(W,idx_global,idx_global);
		llt.compute(W+Kinv);
		L = llt.matrixL(); //no need to extract the triangular matrix here
		//cout<<"L"<<endl<<L<<endl;
		fvis=llt.solve(GetVec(deriv,idx_global)+W*fvis);
		//cout<<"fvis"<<endl<<fvis<<endl;
		for(int w=0;w<idx_global.rows();w++)
		{
			f(idx_global(w))=fvis(w);	
		}
				//cout<<"f"<<endl<<f<<endl;
		loglike =log_likelihood(f, sigma, all_pairs, idx_global_1, idx_global_2,M, N);
		psi_new = loglike - 0.5 * fvis.transpose()*Kinv*fvis; 
	}





}




double Gppe::log_likelihood(VectorXd f,double sigma, TypePair all_pairs,VectorXd idx_global_1,VectorXd idx_global_2,int M, int N)
{
	M=all_pairs.rows();
	VectorXd idx_1,idx_2,z;
	double loglike=0;
	for(int j=0;j<M;j++)
		{	
			if(all_pairs(j).rows()==0)
				continue;
			dsp("hello1");
			idx_1=ind2global(all_pairs(j).col(0),j,N);
			idx_2=ind2global(all_pairs(j).col(1),j,N);
			dsp("hello2");
			dsp(f,"f");
			dsp(idx_1,"idx_1");
			dsp(idx_2,"idx_2");
			z=(GetVec(f,idx_1)-GetVec(f,idx_2))/sigma;
			z=normcdf(z);
			loglike+= log(z.array()).sum();
			dsp("hello3");
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
	
	MatrixXd Deriv2=MatrixXd::Zero(n,n);
	
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

	Deriv2=SetMatGenIdx(Deriv2,ind,-val);

	ind_trans=sub2ind(n,n,idx_global_2,idx_global_1);

	Deriv2=SetMatGenIdx(Deriv2,ind_trans,-val);

	Deriv2=SetMatGenIdx(Deriv2,all_diag_idx,GetMatGenIdx(Deriv2,all_diag_idx)+Get_Cumulative_Val(idx_global_1, val, n));

	Deriv2=SetMatGenIdx(Deriv2,all_diag_idx,GetMatGenIdx(Deriv2,all_diag_idx)+Get_Cumulative_Val(idx_global_2, val, n));

	return Deriv2;
}
