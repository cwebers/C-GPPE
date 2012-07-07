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


//const double CLearner::negative_marginal_log_likelihood(const column_vector & dltheta)
//{
//    VectorXd theta = DlibtoEigen(dltheta);
//    VectorXd theta_x, theta_t;
//    double sigma;
//    GetTheta(theta_x, theta_t, sigma, theta);
//    covt->SetTheta(theta_t);
//    covx->SetTheta(theta_x);
//    CGppe g = CGppe(covt, covx);
//    g.Approx_CGppe_Laplace( theta_x, theta_t, sigma,
//                           t, x, train_pairs, idx_global, idx_global_1, idx_global_2, ind_t, ind_x, M, N);
//
//    double cond_loglike = g.log_likelihood(sigma, train_pairs, idx_global_1, idx_global_2, M, N);
//    VectorXd fvis = GetVec(g.Getf(), idx_global);
//    MatrixXd L = g.Getllt().matrixU();
//    double margl = -0.5 * (-log(g.GetKinv().determinant()) + 2 * (log(g.GetL().diagonal().array()).sum()))
//                   - 0.5 * fvis.transpose() * g.GetKinv() * fvis + cond_loglike;
//
//    cout << "-margl :" << -margl << endl;
//
//    return -margl;
//
//
//}

// Calculate function at point
const double CLearner::negative_marginal_log_likelihood(const column_vector & point)
{
    VectorXd theta = DlibtoEigen(point);
    VectorXd theta_x, theta_t;
    double sigma;
    GetTheta(theta_x, theta_t, sigma, theta);
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
    dsp(-margl, "nl");
    return -margl;

    //  return this->negative_marginal_log_likelihood(arg);
}


//VectorXd CLearner::gradient_negative_marginal_loglikelihood(VectorXd theta)
column_vector gradient_negative_marginal_loglikelihood(const column_vector & point){
    
    VectorXd theta = DlibtoEigen(point);
    
    return EigentoDlib(theta);
}

//{
//    VectorXd theta_x, theta_t;
//    double sigma;
//    VectorXd dtheta_kt, dtheta_kx;
//    GetTheta(theta_x, theta_t, sigma, theta);
//    covt->SetTheta(theta_t);
//    covx->SetTheta(theta_x);
//    CGppe g = CGppe(covt, covx);
//    VectorXd all_diag_idx;
//    int n = M * N;
//    all_diag_idx = sub2ind(n, n, Nfirst(n), Nfirst(n));
//
//    dtheta_kt = VectorXd::Zero(theta_t.rows());
//    dtheta_kx = VectorXd::Zero(theta_x.rows());
//
//
//    // Laplace approximation to the posterior p(f | D)
//    g.Approx_CGppe_Laplace( theta_x, theta_t, sigma,
//                           t, x, train_pairs, idx_global, idx_global_1, idx_global_2, ind_t, ind_x, M, N);
//
//    VectorXd deriv_loglike_vis = g.deriv_log_likelihood_CGppe_fast( sigma, train_pairs, idx_global_1, idx_global_2, M, N);
//    deriv_loglike_vis = GetVec(deriv_loglike_vis, idx_global);
//    MatrixXd Kt = covt->ComputeGrandMatrix(t);
//
//    //Let's Compute things
//    MatrixXd Cinv = (g.Getllt()).solve(g.GetKinv());
//    VectorXd fvis = GetVec(g.Getf(), idx_global);
//    VectorXd alpha = fvis.transpose() * g.GetKinv();
//
//    // We compute the explicit derivative theta_t and theta_x
//
//    double L_thetat = theta_t.rows();
//
//}

/*
VectorXd Learn::Optimize(VectorXd theta_first)
{


        // Now lets try doing it again with a different starting point and the version
        // of find_min() that doesn't require you to supply a derivative function.
        // This version will compute a numerical approximation of the derivative since
        // we didn't supply one to it.
        column_vector starting_point;
        starting_point = EigentoDlib(theta_first);
        find_min_using_approximate_derivatives(bfgs_search_strategy(),
                                               objective_delta_stop_strategy(1e-7),
                                               &negative_marginal_log_likelihood, starting_point, -1);
        // Again the correct minimum point is found and stored in starting_point
        cout << starting_point << endl;

}
*/



////  Calculate function at point
//double CLearner::operator() ( const column_vector& point) const
//{
//    VectorXd theta = DlibtoEigen(point);
//    VectorXd theta_x, theta_t;
//    double sigma;
//    GetTheta(theta_x, theta_t, sigma, theta);
//    sigma = exp(sigma);
//    covt->SetTheta(theta_t);
//    covx->SetTheta(theta_x);
//    CGppe g = CGppe(covt, covx);
//    g.Approx_CGppe_Laplace( theta_x, theta_t, sigma,
//                           t, x, train_pairs, idx_global, idx_global_1, idx_global_2, ind_t, ind_x, M, N);
//
//    double cond_loglike = g.log_likelihood(sigma, train_pairs, idx_global_1, idx_global_2, M, N);
//    VectorXd fvis = GetVec(g.Getf(), idx_global);
//    double margl = -0.5 * (-log(g.GetKinv().determinant()) + 2 * (log(g.GetL().diagonal().array()).sum()))
//                   - 0.5 * fvis.transpose() * g.GetKinv() * fvis + cond_loglike;
//    dsp(-margl, "nl");
//    return -margl;
//
//    //  return this->negative_marginal_log_likelihood(arg);
//}




//  Calculate gradient at point
//int CLearner::operator() ( const column_vector & point) const
//{
//    //column_vector result[3];
//    
//    int result = 1;
//    return result;
//}












