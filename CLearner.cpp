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
#include "CLearner.h"

// Constructor
CLearner::CLearner()
{

}


CLearner::CLearner(Covfunc *Covt, Covfunc *Covx,
             MatrixXd T, MatrixXd X, TypePair All_pairs, VectorXd Idx_global, VectorXd Idx_global_1,
             VectorXd Idx_global_2, VectorXd Ind_t, VectorXd Ind_x, int  m, int n)
{
    covt = Covt;
    covx = Covx;
    t = T;
    x = X;
    train_pairs = All_pairs;
    idx_global = Idx_global;
    idx_global_1 = Idx_global_1;
    idx_global_2 = Idx_global_2;
    ind_t = Ind_t;
    ind_x = Ind_x;
    M = m;
    N = n;
}

// Destructor
CLearner::~CLearner()
{
}


// Calculate function at point
const double CLearner::negative_marginal_log_likelihood(const column_vector & point)
{
    VectorXd theta = DlibtoEigen(point);
    VectorXd theta_x, theta_t;
    double sigma;
    GetTheta(theta_t, theta_x, sigma, theta,t.cols(), x.cols());
    sigma = exp(sigma);
    covt->SetTheta(theta_t);
    covx->SetTheta(theta_x);
    CGppe g = CGppe(covt, covx);
    g.Approx_CGppe_Laplace( theta_x, theta_t, sigma,
                           t, x, train_pairs, idx_global, idx_global_1, idx_global_2, ind_t, ind_x, M, N);

    double cond_loglike = g.log_likelihood(sigma, train_pairs, idx_global_1, idx_global_2, M, N);
    VectorXd fvis = GetVec(g.Getf(), idx_global);
    double margl = -0.5 * (-log(g.GetKinv().determinant()) + 2 * (log(g.GetL().diagonal().array()).sum()))
                   - 0.5 * fvis.transpose() * g.GetKinv() * fvis + cond_loglike;
    //dsp(-margl, "nl");
    return -margl;

}


column_vector CLearner::gradient_negative_marginal_loglikelihood(const column_vector & point)
{    
    VectorXd theta = DlibtoEigen(point);
	MatrixXd dWdf;
    VectorXd theta_x, theta_t;
    double sigma;
    MatrixXd pairs;
    VectorXd dtheta_kt, dtheta_kx;
    GetTheta(theta_t, theta_x, sigma, theta, t.cols(), x.cols());
    sigma=exp(sigma);
    covt->SetTheta(theta_t);
    covx->SetTheta(theta_x);
    CGppe g = CGppe(covt, covx);
    VectorXd all_diag_idx;
    int n = M * N;
    all_diag_idx = sub2ind(n, n, Nfirst(n), Nfirst(n));
	
    dtheta_kt = VectorXd::Zero(theta_t.rows());
    dtheta_kx = VectorXd::Zero(theta_x.rows());

    // Laplace approximation to the posterior p(f | D)
    g.Approx_CGppe_Laplace( theta_x, theta_t, sigma,
                           t, x, train_pairs, idx_global, idx_global_1, idx_global_2, ind_t, ind_x, M, N);

    VectorXd deriv_loglike_vis = g.deriv_log_likelihood_CGppe_fast( sigma, train_pairs, idx_global_1, idx_global_2, M, N);
    deriv_loglike_vis = GetVec(deriv_loglike_vis, idx_global);
    MatrixXd Kt = covt->ComputeGrandMatrix(t);
    
    
    
    //Let's Compute things
    MatrixXd Cinv = (g.Getllt()).solve(g.GetKinv());
    VectorXd fvis = GetVec(g.Getf(), idx_global);
    VectorXd alpha = fvis.transpose() * g.GetKinv();
    
    // We compute the explicit derivative theta_t and theta_x

    double L_thetat = theta_t.rows();
    VecMat dK_dthetat(L_thetat);
    MatrixXd dKt_dthetat;
	    double inter;

    	for(int i=0;i<L_thetat;i++)
    	{
    		dKt_dthetat=covt->Computegrad(t,i);
    		dK_dthetat(i)=GetMat(dKt_dthetat,ind_t,ind_t).array()*GetMat(g.GetKx(),ind_x,ind_x).array();
    		dtheta_kt(i)=-0.5*(Cinv*dK_dthetat(i)*g.GetW()).trace();
    		inter=alpha.transpose()*dK_dthetat(i)*alpha;
    		inter=0.5*inter;
    		dtheta_kt(i)+=inter;
    	}

    double L_thetax = theta_x.rows();
    VecMat dK_dthetax(L_thetax);
    MatrixXd dKx_dthetax;
    	for(int i=0;i<L_thetax;i++)
    	{
    		dKx_dthetax=covx->Computegrad(x,i);
    		dK_dthetax(i)=GetMat(Kt,ind_t,ind_t).array()*GetMat(dKx_dthetax,ind_x,ind_x).array();
    		dtheta_kx(i)=-0.5*(Cinv*dK_dthetax(i)*g.GetW()).trace();
    		inter=alpha.transpose()*dK_dthetax(i)*alpha;
    		inter=0.5*inter;
    		dtheta_kx(i)+=inter;
    	}


	// explicite derivatives wrt sima
	// Need to compute dW_sigma (which is actually dW/dtheta_sigma)
	MatrixXd dWdsigma;
	double dloglike_dsigma, dtheta_sigma;
	get_dsigma(dWdsigma, dloglike_dsigma, g.Getf(), sigma, train_pairs, M, N);
	//dsp(dWdsigma,"dWdsigma");
	//dsp(dloglike_dsigma,"dloglike_dsigma");
	dWdsigma=GetMat(dWdsigma,idx_global, idx_global);
	dtheta_sigma= -0.5*((g.Getllt()).solve(dWdsigma)).trace()+ dloglike_dsigma;
	//dsp(dtheta_sigma,"dtheta_sigma");
	//Compute implicit derivatives here
				

	VectorXd dlogp_dsigma= get_dlogp_dsigma(dWdsigma, dloglike_dsigma, g.Getf(), sigma, train_pairs, M, N);
	VectorXd dfdsigma_vis=(g.Getllt()).solve(GetVec(dlogp_dsigma, idx_global)); //dfdsigma at "observed" values
	//dsp(dfdsigma_vis,"dfdsigma_vis");

	double dtheta_sigma_implicit=0;
	int ptr_ind_t, ptr_ind_x;
	
	for(int i=0;i<idx_global.rows();i++)
	{
		ptr_ind_t=ind_t(i); //indices of current f_{o}
		ptr_ind_x=ind_x(i);
		pairs=train_pairs(ptr_ind_t);

		dWdf=get_dWdf(all_diag_idx, g.Getf(), ptr_ind_t, ptr_ind_x, sigma, pairs, M, N);
				//dsp(dWdf,"dWdf");
		dWdf=GetMat(dWdf,idx_global, idx_global);

		double tmp_val= -0.5*(g.Getllt().solve(dWdf)).trace();

		double val=tmp_val*dfdsigma_vis(i); //tmp_val*dF_dsigma_{i} as dfdsigma only contains "observed values"
		dtheta_sigma_implicit+=val;
		VectorXd df_dtheta;

		for(int k=0;k<theta_t.rows();k++)
		{
			df_dtheta=Cinv*dK_dthetat(k)*deriv_loglike_vis;
			dtheta_kt(k)+=tmp_val*df_dtheta(i);
		}

 		for(int k=0;k<theta_x.rows();k++)
		{
			df_dtheta=Cinv*dK_dthetax(k)*deriv_loglike_vis;
			dtheta_kx(k)+=tmp_val*df_dtheta(i);
		}
	}

	dtheta_sigma+=dtheta_sigma_implicit;
	
	VectorXd grad_theta=-1*concatTheta( dtheta_kt, dtheta_kx, dtheta_sigma);
	//dsp(grad_theta,"grad_theta");
	return EigentoDlib(grad_theta);

}














