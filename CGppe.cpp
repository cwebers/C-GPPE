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
#include "CGppe.h"

CGppe::CGppe(Covfunc *Covfunc_t, Covfunc *Covfunc_x)
{
    covfunc_t = Covfunc_t;
    covfunc_x = Covfunc_x;
}



CGppe::CGppe(const CGppe & g)
{
    f = g.f;
    Kx = g.Kx;
    Kinv = g.Kinv;
    W = g.W;
    L = g.L;
}

CGppe::CGppe(VectorXd fnew, MatrixXd Kxnew, MatrixXd Kinvnew, MatrixXd Wnew, MatrixXd Lnew)
{
    f = fnew;
    Kx = Kxnew;
    Kinv = Kinvnew;
    W = Wnew;
    L = Lnew;

}

CGppe::~CGppe()
{

}


VectorXd CGppe::Getf()
{
    return f;
}

MatrixXd CGppe::GetW()
{
    return W;
}

MatrixXd CGppe::GetL()
{
    return L;
}

LLT<MatrixXd> CGppe::Getllt()
{
    return llt;
}

MatrixXd CGppe::GetKinv()
{
    return Kinv;
}

MatrixXd CGppe::GetKx()
{
    return Kx;
}

VectorXd CGppe::Getmustar()
{
    return mustar;
}

VectorXd CGppe::Getvarstar()
{
    return varstar;
}

double CGppe::Getp()
{
    return p;
}


double CGppe::get_fbest(int N)
{
    VectorXd ftest;
    double fbest;

    ftest = f.segment(f.rows() - N, N - 1);
    fbest = ftest.maxCoeff();
    if (fbest != fbest)
        fbest = f.maxCoeff();
    return fbest;
}

double CGppe::maximum_expected_improvement(const VectorXd & theta_t, const VectorXd& theta_x, const double& sigma,
        const MatrixXd& t, const MatrixXd & x, const VectorXd& idx_global, const VectorXd& ind_t, const VectorXd& ind_x, MatrixXd tstar, int N, double fbest)
{
    VectorXd idx_xstar=Nfirst(N);
    int Kt_ss = 1;
    double  mei;
    MatrixXd Kx_star, Kx_star_star, kstar, Kss, Css;
    MatrixXd Kt_star = covfunc_t->Compute(t, tstar);
	//dsp(GetKinv(),"Kinv");


    Kx_star = GetMatRow(Kx, idx_xstar.transpose()); //maybe need some transpose?

    Kx_star_star = GetMat(Kx, idx_xstar.transpose(), idx_xstar.transpose()); // test to test
    kstar = Kron(Kt_star, Kx_star);

    kstar = GetMatRow(kstar, idx_global);
    Kss = Kt_ss * Kx_star_star;


    mustar = kstar.transpose() * Kinv * GetVec(f, idx_global);
    Css    = Kss - kstar.transpose() * W * llt.solve(Kinv * kstar);
    varstar = Css.diagonal();


    VectorXd sigmastar = sqrt(varstar.array());
    VectorXd z = (fbest - mustar.array()) / sigmastar.array();
    VectorXd pdfval = normpdf(z);
    VectorXd cdfval = normcdf(z);
    VectorXd inter = z.array() * (1 - cdfval.array());
    VectorXd el = sigmastar.cwiseProduct(inter - pdfval);

	el=-1*el;
    mei = el.maxCoeff();
    //dsp(mei,"mei");
    return mei;
}

double CGppe::expected_voi(const VectorXd & theta_x, const VectorXd& theta_t, const double& sigma,
                          const MatrixXd& t, const MatrixXd & x, TypePair train_pairs, VectorXd& idx_global, VectorXd& ind_t, VectorXd& ind_x, MatrixXd test_pair, double fbest, double p_12)
{

    int M = t.rows();
    int N = x.rows();
    CGppe gnew= CGppe(new CovSEard(),new CovSEard());
    VectorXd idx_global_1, idx_global_2;
    MatrixXd tstar = t.row(M-1);
    double p_21, mei_12, mei_21;
    p_21 = 1. - p_12;
	train_pairs(M-1)=MatAdd(train_pairs(M-1),test_pair);
    compute_global_index(idx_global_1, idx_global_2, train_pairs, N);
    unique(idx_global, idx_global_1, idx_global_2);
    ind2sub(ind_x, ind_t, N, M, idx_global);
    gnew.Approx_CGppe_Laplace( theta_x, theta_t, sigma,
                         t, x, train_pairs, idx_global, idx_global_1, idx_global_2, ind_t, ind_x, M, N);

    mei_12 = gnew.maximum_expected_improvement(theta_t, theta_x, sigma, t, x, idx_global, ind_t, ind_x, tstar, N, fbest);

	//dsp(mei_12,"mei12");
    //recomputation
    fliplr(test_pair);
	train_pairs(M-1).bottomRows(1)=test_pair;
    compute_global_index(idx_global_1, idx_global_2, train_pairs, N);
    unique(idx_global, idx_global_1, idx_global_2);
    ind2sub(ind_x, ind_t, N, M, idx_global);


    gnew.Approx_CGppe_Laplace( theta_x, theta_t, sigma,
                         t, x, train_pairs, idx_global, idx_global_1, idx_global_2, ind_t, ind_x, M, N);

    mei_21 = gnew.maximum_expected_improvement(theta_t, theta_x, sigma, t, x, idx_global, ind_t, ind_x, tstar, N, fbest);

    return (p_12*mei_12 + p_21*mei_21) ;
}

void CGppe::Elicit( const VectorXd & theta_x, const VectorXd& theta_t, const double& sigma, const MatrixXd& train_t, const MatrixXd &x, TypePair & train_pairs
                   , const MatrixXd & test_t, int test_user_idx, MatrixXd  idx_pairs, int  Maxiter, const  TypePair& Oracle , MatrixXd& F)
{
    train_pairs.conservativeResize(train_pairs.rows()+1);
    int N = x.rows();
    int Mtrain = train_t.rows();
    int M = Mtrain + 1;
    int Npairs = idx_pairs.rows();
    int Lgood;
    VectorXd vrand, idx_good;
    //VectorXd is_selected(Npairs);
    Matrix<bool, Dynamic, 1> is_selected(Npairs);
    is_selected.fill(false);
    VectorXd loss = VectorXd::Zero(Maxiter + 1);
    double loss_current;
    VectorXd evoi(Npairs), ind_t, ind_x;
    VectorXd idx_global_1, idx_global_2, idx_global;
    compute_global_index(idx_global_1, idx_global_2, train_pairs, N);
    unique(idx_global, idx_global_1, idx_global_2);
    ind2sub(ind_x, ind_t, N, M, idx_global);
    bool stop = false;
    double foo, val;
    int count = 0;
    MatrixXd new_pair;
    MatrixXd t;
	VectorXd p_12(Npairs);
    t.resize(M, train_t.cols());
    t << train_t, test_t;
    for (int iter = 0;iter <= Maxiter;iter++)
    {
        Approx_CGppe_Laplace( theta_x, theta_t, sigma,
                             t, x, train_pairs, idx_global, idx_global_1, idx_global_2, ind_t, ind_x, Mtrain, N);

        Predictive_Utility_Distribution(t, test_t, N, idx_global );
		//dsp(mustar,"mustar");
		//dsp(varstar,"varstar");
        std::ptrdiff_t best_item_idx;
        foo = mustar.maxCoeff(&best_item_idx);
        double fbest = get_fbest(N);
        //dsp(fbest,"fbest");
        MatrixXd test_pair;
        for (int i = 0;i < Npairs;i++)
        {
        	Predict_CGppe_Laplace(sigma, t, x,  idx_global, ind_t, ind_x, t.row(M-1), idx_pairs.row(i));
        	p_12(i)=p;
        }
        
        for (int i = 0;i < Npairs;i++)
        {
            if (is_selected(i))
            {
                evoi(i) = INT_MIN;
                continue;
            }

            test_pair = idx_pairs.row(i);
			evoi(i) = expected_voi(theta_x, theta_t, sigma, t, x, train_pairs, idx_global, ind_t, ind_x, test_pair, fbest,p_12(i));
        }
					//dsp(evoi,"evoi");

        std::ptrdiff_t query_idx;
        val = evoi.maxCoeff(&query_idx);
        idx_good = find(evoi, val);
        //dsp(val,"val");
         Lgood = idx_good.rows(); 
         
    	if ( Lgood > 1) 
    	{
        	vrand = randperm(Lgood);
        	cout<<"Solving clashes at random"<<endl;
        	query_idx = idx_good(vrand(0));
        }
        
        is_selected(query_idx) = true;
		//dsp(query_idx,"queryidx");
        new_pair = make_query_toydata(Oracle, query_idx, test_user_idx);
        //adding the new pair
        train_pairs(M-1)=MatAdd(train_pairs(M-1),new_pair);
        compute_global_index(idx_global_1, idx_global_2, train_pairs, N);

        unique(idx_global, idx_global_1, idx_global_2);

        ind2sub(ind_x, ind_t, N, M, idx_global);

    
     //Computes the loss of making a recommendation at this point
     loss_query_toydata(loss_current, F, stop, test_user_idx, best_item_idx);
     loss(iter)=loss_current;

        count++;
        cout << "Query " << count << "[" << new_pair(0) << " " << new_pair(1) << "] done, Recommended Item= " << best_item_idx << ", loss=" << loss(iter) << endl;
    }
}


void CGppe::Make_Predictions_New_User(const VectorXd & theta_x, const VectorXd& theta_t, double& sigma, const MatrixXd& train_t, const MatrixXd &x, const TypePair & train_pairs,
                                     const VectorXd & idx_global, const VectorXd& idx_global_1, const VectorXd& idx_global_2,
                                     const VectorXd& ind_t, const VectorXd& ind_x, const MatrixXd & test_t, const MatrixXd& idx_pairs, const VectorXd& ftrue, const VectorXd& ytrue)
{
    int N = x.rows();
    int Mtrain = train_t.rows();
    int Npairs = idx_pairs.rows();
    VectorXd fstar;
    MatrixXd pair;
    VectorXd P = VectorXd::Zero(Npairs);
    VectorXd ypred = VectorXd::Zero(Npairs);
    VectorXd sum = VectorXd::Zero(N);
    VectorXd count = VectorXd::Zero(N);

    Approx_CGppe_Laplace( theta_x, theta_t, sigma,
                         train_t, x, train_pairs, idx_global, idx_global_1, idx_global_2, ind_t, ind_x, Mtrain, N);

    for (int i = 0;i < Npairs;i++)
    {
        pair = idx_pairs.row(i);


        Predict_CGppe_Laplace(sigma, train_t, x, idx_global, ind_t, ind_x,
                             test_t, pair);
        P(i) = p;
        sum(pair(0)) += mustar(0);
        count(pair(0)) += 1;
        sum(pair(1)) += mustar(1);
        count(pair(1)) += 1;

    }

    for (int i = 0;i < P.rows();i++)
    {
        if (P(i) > 0.5)
            ypred(i) = 1;
        else
        	ypred(i)=0;
    }

    fstar = sum.array() / count.array();
	dsp(fstar,"fstar");
    cout << endl << endl << "error =  " << (GetDiff(ytrue, ypred)).sum() / ytrue.rows()<<endl;
    // need for a plot function here ?
}





void CGppe::Predictive_Utility_Distribution(MatrixXd t, MatrixXd tstar, int N, VectorXd idx_global)
{
    int Kt_ss = 1;
    VectorXd idx_xstar(N);
    MatrixXd Kstar, Kx_star_star, Kx_star, Kss, Css, Kt_star;
    for (int i = 0;i < N;i++)
    {
        idx_xstar(i) = i;
    }

    Kt_star = covfunc_t->Compute(t, tstar);
    Kx_star = GetMatRow(Kx, idx_xstar);//need to check for tranpose later?
    Kx_star_star = GetMat(Kx, idx_xstar, idx_xstar);

    Kstar = Kron(Kt_star, Kx_star);

    Kstar = GetMatRow(Kstar, idx_global);

    Kss = Kt_ss * Kx_star_star;
    mustar = Kstar.transpose() * Kinv * GetVec(f, idx_global);
    Css = Kss - Kstar.transpose() * W * llt.solve(Kinv * Kstar);
    varstar = Css.diagonal();
}


void CGppe::Predict_CGppe_Laplace(double sigma, MatrixXd t, MatrixXd x, VectorXd idx_global, VectorXd ind_t, VectorXd ind_x,
                                MatrixXd tstar, MatrixXd test_pair)
{

    int Kt_ss = 1;
    double sigma_star, val;
    MatrixXd Kx_star, Kx_star_star, kstar, Kss, Css;
    MatrixXd Kt_star = covfunc_t->Compute(t, tstar);

    Kx_star = GetMatRow(Kx, test_pair.transpose()).transpose(); //maybe need some transpose?

    Kx_star_star = GetMat(Kx, test_pair.transpose(), test_pair.transpose()); // test to test
    kstar = Kron(Kt_star, Kx_star);
    kstar = GetMatRow(kstar, idx_global);

    Kss = Kt_ss * Kx_star_star;

    mustar = kstar.transpose() * Kinv * GetVec(f, idx_global);
    Css    = Kss - kstar.transpose() * W * llt.solve(Kinv * kstar);

    sigma_star = sqrt(Css(0, 0) + Css(1, 1) - 2 * Css(0, 1) + pow(sigma, 2));
    val = ( mustar(0) - mustar(1) ) / sigma_star;
    p   = normcdf(val);

}




void CGppe::Approx_CGppe_Laplace(const VectorXd & theta_x, const VectorXd& theta_t, const double& sigma, const MatrixXd& t, const MatrixXd &x, const TypePair & all_pairs,
                               const VectorXd & idx_global, const VectorXd& idx_global_1, const VectorXd& idx_global_2,
                               const VectorXd& ind_t, const VectorXd& ind_x, int M, int N)
{
    //Parameters function initialization
    double eps = 1E-6, psi_new, psi_old;
    M = all_pairs.rows();
    int n = M * N;
    f = VectorXd::Zero(n);
    VectorXd fvis = VectorXd::Zero(idx_global.rows());
    VectorXd deriv;
    double loglike = 0;

    covfunc_t->SetTheta(theta_t);
    covfunc_x->SetTheta(theta_x);

    MatrixXd Kt = covfunc_t->ComputeGrandMatrix(t);
    Kx = covfunc_x->ComputeGrandMatrix(x);

    MatrixXd K = GetMat(Kt, ind_t, ind_t).array() * GetMat(Kx, ind_x, ind_x).array();

    loglike = log_likelihood( sigma, all_pairs, idx_global_1, idx_global_2, M, N);
    Kinv = K.inverse();
    psi_new = loglike - 0.5 * fvis.transpose() * Kinv * fvis;
    psi_old = INT_MIN;
    while ((psi_new - psi_old) > eps)
    {
        psi_old = psi_new;
        deriv = deriv_log_likelihood_CGppe_fast( sigma, all_pairs, idx_global_1, idx_global_2, M, N);
        W = -deriv2_log_likelihood_CGppe_fast(sigma, all_pairs, idx_global_1, idx_global_2, M, N);
        W = GetMat(W, idx_global, idx_global);
        llt.compute(W + Kinv);
        L = llt.matrixL(); //no need to extract the triangular matrix here
        fvis = llt.solve(GetVec(deriv, idx_global) + W * fvis);
        for (int w = 0;w < idx_global.rows();w++)
        {
            f(idx_global(w)) = fvis(w);
        }
        loglike = log_likelihood( sigma, all_pairs, idx_global_1, idx_global_2, M, N);
        psi_new = loglike - 0.5 * fvis.transpose() * Kinv * fvis;
    }





}




double CGppe::log_likelihood(double sigma, TypePair all_pairs, VectorXd idx_global_1, VectorXd idx_global_2, int M, int N)
{
    M = all_pairs.rows();
    VectorXd idx_1, idx_2, z;
    double loglike = 0;
    for (int j = 0;j < M;j++)
    {
        if (all_pairs(j).rows() == 0)
            continue;
        idx_1 = ind2global(all_pairs(j).col(0), j, N);
        idx_2 = ind2global(all_pairs(j).col(1), j, N);
        z = (GetVec(f, idx_1) - GetVec(f, idx_2)) / sigma;
        z = normcdf(z);
        loglike += log(z.array()).sum();
    }
    return loglike;
}

VectorXd CGppe::deriv_log_likelihood_CGppe_fast(double sigma, const  TypePair& all_pairs, VectorXd idx_global_1, VectorXd idx_global_2, int M, int N)
{
    VectorXd deriv_loglike, z, cdf_val, pdf_val, val;
    // test if idx vectors are empty ?
    M = all_pairs.rows();
    int n = M * N;
    z = (GetVec(f, idx_global_1) - GetVec(f, idx_global_2)) / sigma;
    cdf_val = normcdf(z);
    pdf_val = normpdf(z);
    val = pow(sigma,-1) * (pdf_val.cwiseQuotient(cdf_val));
    return Get_Cumulative_Val(idx_global_1, val, n) - Get_Cumulative_Val(idx_global_2, val, n);
}

MatrixXd CGppe::deriv2_log_likelihood_CGppe_fast(double sigma, TypePair all_pairs, VectorXd idx_global_1, VectorXd idx_global_2, int M, int N)
{
    VectorXd deriv_loglike, z, cdf_val, pdf_val, val, ratio, all_diag_idx, ind, ind_trans;

    M = all_pairs.rows();

    int n = M * N;

    VectorXd consec(n);

    MatrixXd Deriv2 = MatrixXd::Zero(n, n);

    for (int i = 0;i < n;i++)
    {
        consec(i) = i;
    }
    all_diag_idx = sub2ind(n, n, consec, consec);

    z = (GetVec(f, idx_global_1) - GetVec(f, idx_global_2)) / sigma;

    cdf_val = normcdf(z);

    pdf_val = normpdf(z);

    ratio = pdf_val.array() / cdf_val.array();

    val = -(1. / pow(sigma, 2)) * (ratio.array() * (z + ratio).array());

    ind = sub2ind(n, n, idx_global_1, idx_global_2);

    Deriv2 = SetMatGenIdx(Deriv2, ind, -val);

    ind_trans = sub2ind(n, n, idx_global_2, idx_global_1);

    Deriv2 = SetMatGenIdx(Deriv2, ind_trans, -val);

    Deriv2 = SetMatGenIdx(Deriv2, all_diag_idx, GetMatGenIdx(Deriv2, all_diag_idx) + Get_Cumulative_Val(idx_global_1, val, n));

    Deriv2 = SetMatGenIdx(Deriv2, all_diag_idx, GetMatGenIdx(Deriv2, all_diag_idx) + Get_Cumulative_Val(idx_global_2, val, n));

    return Deriv2;
}
