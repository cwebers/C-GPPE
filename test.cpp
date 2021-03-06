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

int testmat()
{
    //for measuring running time
    clock_t start, end;
    double elapsed;
    start = clock();
    
    MatrixXd a(5000,5000), b(5000,5000),c;
    a.setRandom();
    b.setRandom();
    c=a*b;
    
	end = clock();
    elapsed = ((double)end - start) / CLOCKS_PER_SEC;
    cout << "Elapsed Time :" << elapsed << endl;

}


int General()
{
	//Here we are testing Prediction, Elicitation and Optimisation in one shot
	

    
    //for measuring running time
   // clock_t start, end;
    double elapsed;
    //start = clock();
	//declaring the data
	int N,M;
	int Maxiter=10;
	TypePair train_pairs, Oracle;
	double logsigma=-2.3025;
	double sigma=0.1;
 	VectorXd idx_global_1, idx_global_2, idx_global, ind_t, ind_x;
 	MatrixXd pairs, test_t, t, x,idx_pairs, F;
 	VectorXd theta_x;
 	VectorXd theta_t;
 	VectorXd theta;
 	VectorXd ftrue, ytrue;
 	MatrixXd test_pairs;
 	MatrixXd train_t;
 	CGppe g=CGppe(new CovSEard(),new CovSEard());

 	//assigning the Data
 	Generate(idx_pairs, t, x, Oracle, train_pairs, F, new CovSEard(), new CovSEard(), theta_t, theta_x, M, N, ftrue,ytrue, test_pairs, test_t 
 			,train_t);
 	//Computing the indexes
 	compute_global_index(idx_global_1, idx_global_2, train_pairs, N);
    unique(idx_global, idx_global_1, idx_global_2);
    ind2sub(ind_x, ind_t, N, M, idx_global);
     int test_user_idx=M-1;
    theta=concatTheta(theta_t, theta_x,logsigma);
    dsp(M,"Number of Users");
    dsp(N,"Number of Items");
    dsp(theta_t.rows()-1,"Dimension of Users");
    dsp(theta_x.rows()-1,"Dimension of Items");
    int Mtrain=M-1;
    
	     g.Make_Predictions_New_User(theta_x, theta_t, sigma, t, x, train_pairs,
                                   idx_global, idx_global_1, idx_global_2,
                                   ind_t, ind_x, test_t, idx_pairs, ftrue, ytrue);
                                   
	g.Elicit(theta_x, theta_t, sigma,  train_t, x, train_pairs, test_t, 
		test_user_idx, idx_pairs, Maxiter,  Oracle, F);
	
	
	
	/*
   // VectorXd theta_first = theta;
    VectorXd theta_first = VectorXd::Zero(theta.rows());
    theta_first(theta.rows()-1)=logsigma;

    CLearner learner = CLearner(new CovSEard(), new CovSEard(),
                          train_t, x, train_pairs, idx_global, idx_global_1, idx_global_2, ind_t, ind_x, Mtrain, N);





    // Now lets try doing it again with a different starting point and the version
    // of find_min() that doesn't require you to supply a derivative function.
    // This version will compute a numerical approximation of the derivative since
    // we didn't supply one to it.
    column_vector starting_point;
    starting_point = EigentoDlib(theta_first);

        cout << "Difference between analytic derivative and numerical approximation of derivative: " 
              << length(derivative(CNLL_Function(learner))(starting_point) - learner.gradient_negative_marginal_loglikelihood(starting_point)) << endl;


// find_min_using_approximate_derivatives(bfgs_search_strategy(),
// objective_delta_stop_strategy(1e-7),
// learner,
// starting_point, INT_MIN);
 
    find_min(bfgs_search_strategy(),
                                           objective_delta_stop_strategy(1e-7),
                                           CNLL_Function(learner),
                                           CGradNLL_Function(learner),
                                           starting_point, -10);
    
    // Again the correct minimum point is found and stored in starting_point
    cout<<"grad_theta"<<endl;
    cout << starting_point << endl;
	*/
	//end = clock();
    //elapsed = ((double)end - start) / CLOCKS_PER_SEC;
    //cout << "Elapsed Time :" << elapsed << endl;
	return 0;	
}


int Optimisation()
{
		//for measuring running time
    clock_t start, end;
    double elapsed;
    start = clock();
	//declaring the data
	int N,M;
	int Maxiter=10;
	TypePair train_pairs, Oracle;
	double logsigma=-2.3025;
 	VectorXd idx_global_1, idx_global_2, idx_global, ind_t, ind_x;
 	MatrixXd pairs, test_t, t, x,idx_pairs, F;
 	VectorXd theta_x;
 	VectorXd theta_t;
 	VectorXd theta;
 	VectorXd ftrue, ytrue;
 	MatrixXd test_pairs;
 	MatrixXd train_t;
 	CGppe g=CGppe(new CovSEard(),new CovSEard());

 	//assigning the Data
 	Generate(idx_pairs, t, x, Oracle, train_pairs, F, new CovSEard(), new CovSEard(), theta_t, theta_x, M, N, ftrue,ytrue, test_pairs, test_t 
 			,train_t);
 	//Computing the indexes
 	compute_global_index(idx_global_1, idx_global_2, train_pairs, N);
    unique(idx_global, idx_global_1, idx_global_2);
    ind2sub(ind_x, ind_t, N, M, idx_global);
     int test_user_idx=M-1;
    theta=concatTheta(theta_t, theta_x,logsigma);
    
    int Mtrain=M-1;
    
   
   
   
   
   
   // VectorXd theta_first = theta;
    VectorXd theta_first = VectorXd::Zero(theta.rows());
    theta_first(theta.rows()-1)=logsigma;

    CLearner learner = CLearner(new CovSEard(), new CovSEard(),
                          train_t, x, train_pairs, idx_global, idx_global_1, idx_global_2, ind_t, ind_x, Mtrain, N);





    // Now lets try doing it again with a different starting point and the version
    // of find_min() that doesn't require you to supply a derivative function.
    // This version will compute a numerical approximation of the derivative since
    // we didn't supply one to it.
    column_vector starting_point;
    starting_point = EigentoDlib(theta_first);

        cout << "Difference between analytic derivative and numerical approximation of derivative: " 
              << length(derivative(CNLL_Function(learner))(starting_point) - learner.gradient_negative_marginal_loglikelihood(starting_point)) << endl;


// find_min_using_approximate_derivatives(bfgs_search_strategy(),
// objective_delta_stop_strategy(1e-7),
// learner,
// starting_point, INT_MIN);
 
    find_min(bfgs_search_strategy(),
                                           objective_delta_stop_strategy(1e-7),
                                           CNLL_Function(learner),
                                           CGradNLL_Function(learner),
                                           starting_point, -10);
    
    // Again the correct minimum point is found and stored in starting_point
    cout<<"grad_theta"<<endl;
    cout << starting_point << endl;

    end = clock();
    elapsed = ((double)end - start) / CLOCKS_PER_SEC;
    cout << "Elapsed Time :" << elapsed << endl;
	return 0;
}



int Posterior()
{
	//for measuring running time
    clock_t start, end;
    double elapsed;
    start = clock();
	//declaring the data
	int N,M;
	int Maxiter=10;
	TypePair train_pairs, Oracle;
	double sigma=0.1;
 	VectorXd idx_global_1, idx_global_2, idx_global, ind_t, ind_x;
 	MatrixXd pairs, test_t, t, x,idx_pairs, F;
 	VectorXd theta_x;
 	VectorXd theta_t;
 	VectorXd theta;
 	VectorXd ftrue, ytrue;
 	MatrixXd test_pairs;
 	MatrixXd train_t;
 	CGppe g=CGppe(new CovSEard(),new CovSEard());

 	//assigning the Data
 	Generate(idx_pairs, t, x, Oracle, train_pairs, F, new CovSEard(), new CovSEard(), theta_t, theta_x, M, N, ftrue,ytrue, test_pairs, test_t 
 			,train_t);
 	//Computing the indexes
 	compute_global_index(idx_global_1, idx_global_2, train_pairs, N);
    unique(idx_global, idx_global_1, idx_global_2);
    ind2sub(ind_x, ind_t, N, M, idx_global);
    g.Approx_CGppe_Laplace( theta_x, theta_t, sigma,
                           t, x, train_pairs, idx_global, idx_global_1, idx_global_2, ind_t, ind_x, M, N);
    dsp(g.GetW(), "W");
    dsp(g.GetL(), "L");
    dsp(g.GetKinv(), "Kinv");
    dsp(g.Getf(), "f");


    end = clock();
    elapsed = ((double)end - start) / CLOCKS_PER_SEC;
    cout << "Elapsed Time :" << elapsed << endl;
	return 0;
}


int teststring()
{
	TypePair all_pairs;
	all_pairs=InputPair("/Users/christopheroustel/Desktop/C-GPPE/all_pairs");
	dspair(all_pairs,"all_pairs");
	dsp(all_pairs.rows(),"avant");
	all_pairs.conservativeResize(all_pairs.rows()+1);
	dsp(all_pairs.rows(),"après");
	dspair(all_pairs,"all_pairs");

	return 0;
}

int Optimisation_without_derivatives()
{
	//for measuring running time
    clock_t start, end;
    double elapsed;
    start = clock();
	//declaring the data
	int N,M;
	int Maxiter=10;
	TypePair train_pairs, Oracle;
	double logsigma=-2.3025;
 	VectorXd idx_global_1, idx_global_2, idx_global, ind_t, ind_x;
 	MatrixXd pairs, test_t, t, x,idx_pairs, F;
 	VectorXd theta_x;
 	VectorXd theta_t;
 	VectorXd theta;
 	VectorXd ftrue, ytrue;
 	MatrixXd test_pairs;
 	MatrixXd train_t;
 	CGppe g=CGppe(new CovSEard(),new CovSEard());

 	//assigning the Data
 	Generate(idx_pairs, t, x, Oracle, train_pairs, F, new CovSEard(), new CovSEard(), theta_t, theta_x, M, N, ftrue,ytrue, test_pairs, test_t 
 			,train_t);
 	//Computing the indexes
 	compute_global_index(idx_global_1, idx_global_2, train_pairs, N);
    unique(idx_global, idx_global_1, idx_global_2);
    ind2sub(ind_x, ind_t, N, M, idx_global);
     int test_user_idx=M-1;
    theta=concatTheta(theta_t, theta_x,logsigma);
    
    int Mtrain=M-1;
    
   
   
   
   
   
   // VectorXd theta_first = theta;
    VectorXd theta_first = VectorXd::Zero(theta.rows());
    theta_first(theta.rows()-1)=logsigma;

    CLearner learner = CLearner(new CovSEard(), new CovSEard(),
                          train_t, x, train_pairs, idx_global, idx_global_1, idx_global_2, ind_t, ind_x, Mtrain, N);





    // Now lets try doing it again with a different starting point and the version
    // of find_min() that doesn't require you to supply a derivative function.
    // This version will compute a numerical approximation of the derivative since
    // we didn't supply one to it.
    column_vector starting_point;
    starting_point = EigentoDlib(theta_first);

        cout << "Difference between analytic derivative and numerical approximation of derivative: " 
              << length(derivative(CNLL_Function(learner))(starting_point) - learner.gradient_negative_marginal_loglikelihood(starting_point)) << endl;


// find_min_using_approximate_derivatives(bfgs_search_strategy(),
// objective_delta_stop_strategy(1e-7),
// learner,
// starting_point, INT_MIN);
 
    find_min_using_approximate_derivatives(bfgs_search_strategy(),
                                           objective_delta_stop_strategy(1e-7),
                                           CNLL_Function(learner),
                                           starting_point, -10);
    
    // Again the correct minimum point is found and stored in starting_point
    cout << starting_point << endl;
    end = clock();
    elapsed = ((double)end - start) / CLOCKS_PER_SEC;
    cout << "Elapsed Time :" << elapsed << endl;
	return 0;
}

int Prediction()
{
	//for measuring running time
    clock_t start, end;
    double elapsed;
    start = clock();
	//declaring the data
	int N,M;
	int Maxiter=10;
	TypePair train_pairs, Oracle;
	double sigma=0.1;
 	VectorXd idx_global_1, idx_global_2, idx_global, ind_t, ind_x;
 	MatrixXd pairs(1, 2), test_t, t, x,idx_pairs, F;
 	VectorXd theta_x;
 	VectorXd theta_t;
 	VectorXd ftrue, ytrue;
 	MatrixXd test_pairs;
 	MatrixXd train_t;
 	CGppe g=CGppe(new CovSEard(),new CovSEard());

 	//assigning the Data
 	Generate(idx_pairs, t, x, Oracle, train_pairs, F, new CovSEard(), new CovSEard(), theta_t, theta_x, M, N, ftrue,ytrue, test_pairs, test_t 
 			,train_t);
 	//Computing the indexes
 	compute_global_index(idx_global_1, idx_global_2, train_pairs, N);
    unique(idx_global, idx_global_1, idx_global_2);
    ind2sub(ind_x, ind_t, N, M, idx_global);
     int test_user_idx=M-1;
     
     
     
     g.Make_Predictions_New_User(theta_x, theta_t, sigma, t, x, train_pairs,
                                   idx_global, idx_global_1, idx_global_2,
                                   ind_t, ind_x, test_t, idx_pairs, ftrue, ytrue);
	
	
		    end = clock();

    elapsed = ((double)end - start) / CLOCKS_PER_SEC;
    cout << "Elapsed Time :" << elapsed << endl;
	return 0;
}
int Elicitation()
{
	//for measuring running time
    clock_t start, end;
    double elapsed;
    start = clock();
	//declaring the data
	int N,M;
	int Maxiter=10;
	TypePair train_pairs, Oracle;
	double sigma=0.1;
 	VectorXd idx_global_1, idx_global_2, idx_global, ind_t, ind_x;
 	MatrixXd pairs(1, 2), test_t, t, x,idx_pairs, F;
 	VectorXd theta_x;
 	VectorXd theta_t;
 	VectorXd ftrue, ytrue;
 	MatrixXd test_pairs;
 	MatrixXd train_t;
 	CGppe g=CGppe(new CovSEard(),new CovSEard());

 	//assigning the Data
 	Generate(idx_pairs, t, x, Oracle, train_pairs, F, new CovSEard(), new CovSEard(), theta_t, theta_x, M, N, ftrue,ytrue, test_pairs, test_t 
 			,train_t);
 	//Computing the indexes
 	compute_global_index(idx_global_1, idx_global_2, train_pairs, N);
    unique(idx_global, idx_global_1, idx_global_2);
    ind2sub(ind_x, ind_t, N, M, idx_global);
     int test_user_idx=M-1;

	
	g.Elicit(theta_x, theta_t, sigma,  train_t, x, train_pairs, test_t, 
		test_user_idx, idx_pairs, Maxiter,  Oracle, F);

	    end = clock();

    elapsed = ((double)end - start) / CLOCKS_PER_SEC;
    cout << "Elapsed Time :" << elapsed << endl;
	return 0;

}

int testgenerate()
{
	//for measuring running time
    clock_t start, end;
    double elapsed;
    start = clock();
	//declaring the data
	int N,M;
	TypePair train_pairs, Oracle;
	double sigma=0.1;
 	VectorXd idx_global_1, idx_global_2, idx_global, ind_t, ind_x;
 	MatrixXd pairs(1, 2), test_t, t, x,idx_pairs, F;
 	VectorXd theta_x;
 	VectorXd theta_t;
 	VectorXd ftrue, ytrue;
 	MatrixXd test_pairs;
 	MatrixXd train_t;
 	CGppe g=CGppe(new CovSEard(),new CovSEard());

 	//assigning the Data
 	Generate(idx_pairs, t, x, Oracle, train_pairs, F, new CovSEard(), new CovSEard(), theta_t, theta_x, M, N, ftrue,ytrue, test_pairs, test_t 
 			,train_t);
 	//Computing the indexes
 	compute_global_index(idx_global_1, idx_global_2, train_pairs, N);
    unique(idx_global, idx_global_1, idx_global_2);
    ind2sub(ind_x, ind_t, N, M, idx_global);

	
	
	    g.Approx_CGppe_Laplace( theta_x, theta_t, sigma,
                           		t, x, train_pairs, idx_global, idx_global_1, idx_global_2, ind_t, ind_x, M, N);
	
	    end = clock();


    elapsed = ((double)end - start) / CLOCKS_PER_SEC;
    cout << "Elapsed Time :" << elapsed << endl;
	return 0;
}



int testaddrows()
{
	VectorXd theta_t(3), theta_x(4), theta;
	double sigma;
	theta_t<<1,2,3;
	theta_x<<4,5,6,7;
	sigma=8;
theta=concatTheta(theta_t,theta_x,sigma);
dsp(theta,"theta");
	GetTheta(theta_t,theta_x,sigma ,theta, 2, 3);
dsp(theta_t,"theta_t");
dsp(theta_x,"theta_x");
dsp(sigma,"sigma");

	return 0;
}


int testreshape()
{
	VectorXd f(6);
	f<<1,2,3,4,5,6;
	dsp(reshape(f,3,2),"reshape");
	return 0;
}


int testelicit()
{
 //generating the data naively

 int M = 3;
 int N = 2;
 int Maxiter=10;
 int test_user_idx=2;
 double sigma = 0.1;
 CGppe g = CGppe(new CovSEard(), new CovSEard());
 TypePair train_pairs(3), Oracle(3);
 VectorXd idx_global_1(2), idx_global_2(2), idx_global(4), ind_t(4), ind_x(4);
 MatrixXd pairs(1, 2), t(2, 2), x(2, 3), test_t(1, 2);
 VectorXd theta_x = VectorXd::Zero(4);
 VectorXd theta_t = VectorXd::Zero(3);
 VectorXd theta = VectorXd::Zero(8);
 theta(7) = -2.3026;
 t(0, 0) = -0.7258;
 t(0, 1) = -1.9623;
 t(1, 0) = -0.3078;
 t(1, 1) = -0.9332;
 x(0, 0) = 2.4582;
 x(0, 1) = -4.0911;
 x(0, 2) = 1.0004;
 x(1, 0) = 6.1426;
 x(1, 1) = -6.3481;
 x(1, 2) = -4.7591;
 MatrixXd train_t(2,2);
 train_t<<   -0.7258 ,  -1.9623,
   -0.3078 ,  -0.9332;
 MatrixXd F(N,M);
 F<<  1.1780  ,  0.8440,   -1.2463,
   -1.1142  , -1.2890,   -1.8551;
 pairs << 0, 1;
 test_t << 0.2501, 1.4168;
 train_pairs(0) = pairs;
 train_pairs(1) = pairs;
 
 Oracle(0) = pairs;
 Oracle(1) = pairs;
 Oracle(2)=pairs;
 
 VectorXd ftrue(2);
 ftrue<<
   -1.2463,
   -1.8551;

MatrixXd idx_pairs(1,2);
idx_pairs<<0,1;

VectorXd ytrue(1);
	ytrue<<1;

 idx_global_1 << 0, 2;
 idx_global_2 << 1, 3;
 idx_global << 0, 1, 2, 3;

 ind_t << 0, 0, 1, 1;
 ind_x << 0, 1, 0, 1;


g.Elicit(theta_x, theta_t, sigma,  train_t, x, train_pairs, test_t, 
		test_user_idx, idx_pairs, Maxiter,  Oracle, F);


	return 0;
}


int testinput()
{
	string const namefile("/Users/christopheroustel/Desktop/C-GPPE/x.txt");
	MatrixXd z=GetData(namefile);
	
	dsp(z,"z");
	return 0;
} 




int testprediction()
{
 //generating the data naively

 int M = 3;
 int N = 2;
 double sigma = 0.1;
 CGppe g = CGppe(new CovSEard(), new CovSEard());
 TypePair all_pairs(2);
 VectorXd idx_global_1(2), idx_global_2(2), idx_global(4), ind_t(4), ind_x(4);
 MatrixXd pairs(1, 2), t(2, 2), x(2, 3), tstar(1, 2);
 VectorXd theta_x = VectorXd::Zero(4);
 VectorXd theta_t = VectorXd::Zero(3);
 VectorXd theta = VectorXd::Zero(8);
 theta(7) = -2.3026;
 t(0, 0) = -0.7258;
 t(0, 1) = -1.9623;
 t(1, 0) = -0.3078;
 t(1, 1) = -0.9332;
 x(0, 0) = 2.4582;
 x(0, 1) = -4.0911;
 x(0, 2) = 1.0004;
 x(1, 0) = 6.1426;
 x(1, 1) = -6.3481;
 x(1, 2) = -4.7591;
 pairs << 0, 1;
 tstar << 0.2501, 1.4168;
 all_pairs(0) = pairs;
 all_pairs(1) = pairs;
 VectorXd ftrue(2);
 ftrue<<
   -1.2463,
   -1.8551;

MatrixXd idx_pairs(1,2);
idx_pairs<<0,1;

VectorXd ytrue(1);
	ytrue<<1;

 idx_global_1 << 0, 2;
 idx_global_2 << 1, 3;
 idx_global << 0, 1, 2, 3;

 ind_t << 0, 0, 1, 1;
 ind_x << 0, 1, 0, 1;

 g.Make_Predictions_New_User(theta_x, theta_t, sigma, t, x, all_pairs,
                                   idx_global, idx_global_1, idx_global_2,
                                   ind_t, ind_x, tstar, idx_pairs, ftrue, ytrue);
}



int testgradcov()
{
	MatrixXd t(2,2);
	VectorXd theta_t=VectorXd::Zero(3);
    t(0, 0) = -0.7258;
    t(0, 1) = -1.9623;
    t(1, 0) = -0.3078;
    t(1, 1) = -0.9332;
    
    CovSEard a=CovSEard(theta_t);
    dsp(a.Computegrad(t,2),"res");
    
	return 0;
}



int testopt()
{
    //generating the data naively

    int M = 3;
    int N = 2;
    double sigma = 0.1;
    CGppe g = CGppe(new CovSEard(), new CovSEard());
    TypePair all_pairs(2);
    VectorXd idx_global_1(2), idx_global_2(2), idx_global(4), ind_t(4), ind_x(4);
    MatrixXd pairs(1, 2), t(2, 2), x(2, 3), tstar(1, 2);
    VectorXd theta_x = VectorXd::Ones(4);
    VectorXd theta_t = VectorXd::Ones(3);
    VectorXd theta = VectorXd::Ones(8);
    theta(7) = -2.30258509;
    VectorXd theta_first = theta;

    t(0, 0) = -0.7258;
    t(0, 1) = -1.9623;
    t(1, 0) = -0.3078;
    t(1, 1) = -0.9332;
    x(0, 0) = 2.4582;
    x(0, 1) = -4.0911;
    x(0, 2) = 1.0004;
    x(1, 0) = 6.1426;
    x(1, 1) = -6.3481;
    x(1, 2) = -4.7591;

    pairs << 0, 1;
    all_pairs(0) = pairs;
    all_pairs(1) = pairs;


    idx_global_1 << 0, 2;
    idx_global_2 << 1, 3;
    idx_global << 0, 1, 2, 3;
    ind_t << 0, 0, 1, 1;
    ind_x << 0, 1, 0, 1;

    CLearner learner = CLearner(new CovSEard(), new CovSEard(),
                          t, x, all_pairs, idx_global, idx_global_1, idx_global_2, ind_t, ind_x, M, N);

    // Now lets try doing it again with a different starting point and the version
    // of find_min() that doesn't require you to supply a derivative function.
    // This version will compute a numerical approximation of the derivative since
    // we didn't supply one to it.
    column_vector starting_point;
    starting_point = EigentoDlib(theta_first);

// find_min_using_approximate_derivatives(bfgs_search_strategy(),
// objective_delta_stop_strategy(1e-7),
// learner,
// starting_point, INT_MIN);
 
    find_min(bfgs_search_strategy(),
                                           objective_delta_stop_strategy(1e-7),
                                           CNLL_Function(learner),CGradNLL_Function(learner),
                                           starting_point, INT_MIN);
    
    // Again the correct minimum point is found and stored in starting_point
    cout << starting_point << endl;

    return 0;
}

int testcovderiv()
{

    VectorXd theta_x = VectorXd::Zero(4);
    MatrixXd x(2, 3);
    x(0, 0) = 2.4582;
    x(0, 1) = -4.0911;
    x(0, 2) = 1.0004;
    x(1, 0) = 6.1426;
    x(1, 1) = -6.3481;
    x(1, 2) = -4.7591;
    MatrixXd t(2, 3);
    t.fill(1);
    CovSEard a = CovSEard(theta_x);
    dsp(a.ComputeGrandMatrix(x), "res");
    return 0;
}

int testgradnl()
{
 //generating the data naively
	//for measuring running time
    clock_t start, end;
    double elapsed;
    start = clock();
	//declaring the data
	int N,M,Mtrain;
	int Maxiter=10;
	TypePair train_pairs, Oracle;
	double logsigma=-2.3026;
 	VectorXd idx_global_1, idx_global_2, idx_global, ind_t, ind_x;
 	MatrixXd pairs, test_t, t, x,idx_pairs, F;
 	VectorXd theta_x;
 	VectorXd theta_t;
 	VectorXd theta;
 	VectorXd ftrue, ytrue;
 	MatrixXd test_pairs;
 	MatrixXd train_t;
 	Mtrain=M-1;
 	//assigning the Data
 	Generate(idx_pairs, t, x, Oracle, train_pairs, F, new CovSEard(), new CovSEard(), theta_t, theta_x, M, N, ftrue,ytrue, test_pairs, test_t 
 			,train_t);
 	//Computing the indexes
 	compute_global_index(idx_global_1, idx_global_2, train_pairs, N);
    unique(idx_global, idx_global_1, idx_global_2);
    ind2sub(ind_x, ind_t, N, M, idx_global);
    	theta=concatTheta(theta_t, theta_x,logsigma);
 	Mtrain=M-1;

 CLearner learner = CLearner(new CovSEard(), new CovSEard(),
	train_t, x, train_pairs, idx_global, idx_global_1, idx_global_2, ind_t, ind_x, Mtrain, N);
dsp(DlibtoEigen(learner.gradient_negative_marginal_loglikelihood(EigentoDlib(theta))),"gradtheta");
 	    end = clock();


    elapsed = ((double)end - start) / CLOCKS_PER_SEC;
    cout << "Elapsed Time :" << elapsed << endl;
	return 0;
}






int testnl()
{
 //generating the data naively

 int M = 3;
 int N = 2;
 double sigma = 0.1;
 CGppe g = CGppe(new CovSEard(), new CovSEard());
 TypePair all_pairs(2);
 VectorXd idx_global_1(2), idx_global_2(2), idx_global(4), ind_t(4), ind_x(4);
 MatrixXd pairs(1, 2), t(2, 2), x(2, 3), tstar(1, 2);
 VectorXd theta_x = VectorXd::Zero(4);
 VectorXd theta_t = VectorXd::Zero(3);
 VectorXd theta = VectorXd::Zero(8);
 theta(7) = -2.3026;
 t(0, 0) = -0.7258;
 t(0, 1) = -1.9623;
 t(1, 0) = -0.3078;
 t(1, 1) = -0.9332;
 x(0, 0) = 2.4582;
 x(0, 1) = -4.0911;
 x(0, 2) = 1.0004;
 x(1, 0) = 6.1426;
 x(1, 1) = -6.3481;
 x(1, 2) = -4.7591;
 pairs << 0, 1;
 tstar << 0.2501, 1.4168;
 all_pairs(0) = pairs;
 all_pairs(1) = pairs;


 idx_global_1 << 0, 2;
 idx_global_2 << 1, 3;
 idx_global << 0, 1, 2, 3;
 ind_t << 0, 0, 1, 1;
 ind_x << 0, 1, 0, 1;
 CLearner learner = CLearner(new CovSEard(), new CovSEard(),
	t, x, all_pairs, idx_global, idx_global_1, idx_global_2, ind_t, ind_x, M, N);

 dsp(learner.negative_marginal_log_likelihood(EigentoDlib(theta)), "nl");
 return 0;
}


int findvalue()
{
    VectorXd theta(7), theta_x(3), theta_t(3);
    double sigma;
    theta << 1, 2, 3, 4, 5, 6, 7;
    GetTheta(theta_x, theta_t, sigma, theta,3,3);
    dsp(theta_x, "theta_x");
    dsp(theta_t, "theta_t");
    dsp(sigma, "sigma");

    dsp(concatTheta(theta_x, theta_t, sigma), "theta");

    return 0;
}

int testgendata()
{
    int M = 3, N = 2;
    TypePair all_pairs(2);
    VectorXd idx_global_1, idx_global_2, idx_global, ind_t(4), ind_x(4);
    MatrixXd pairs(1, 2), t(2, 2), x(2, 3), tstar(1, 2);
    VectorXd theta_x = VectorXd::Zero(4);
    VectorXd theta_t = VectorXd::Zero(3);
    VectorXd f = VectorXd::Zero(6);
    t(0, 0) = -0.7258;
    t(0, 1) = -1.9623;
    t(1, 0) = -0.3078;
    t(1, 1) = -0.9332;
    x(0, 0) = 2.4582;
    x(0, 1) = -4.0911;
    x(0, 2) = 1.0004;
    x(1, 0) = 6.1426;
    x(1, 1) = -6.3481;
    x(1, 2) = -4.7591;
    pairs << 0, 1;
    tstar << 0.2501, 1.4168;
    all_pairs(0) = pairs;
    all_pairs(1) = pairs;
    compute_global_index(idx_global_1, idx_global_2, all_pairs, N);
    dsp(idx_global_1, "idx_global_1");
    dsp(idx_global_2, "idx_global_2");

    unique(idx_global, idx_global_1, idx_global_2);
    ind2sub(ind_x, ind_t, N, M, idx_global);
    dsp(ind_t, "ind_t");
    dsp(ind_x, "ind_x");


    return 0;
}




int testmatrixmultiplication()
{
    clock_t start, end;
    double elapsed;
    start = clock();
    MatrixXd A(5000, 5000), B(5000, 5000), C;
    //A.setRandom();
    //B.setRandom();
    C = A * B;


    end = clock();
    elapsed = ((double)end - start) / CLOCKS_PER_SEC;
    cout << "Elapsed Time :" << elapsed << endl;
    return 0;
}



int testpredictive_utility()
{
    //generating the data naively
    int M = 3;
    int N = 2;
    double sigma = 0.1;
    CGppe g = CGppe(new CovSEard(), new CovSEard());
    TypePair all_pairs(2);
    VectorXd idx_global_1(2), idx_global_2(2), idx_global(4), ind_t(4), ind_x(4);
    MatrixXd pairs(1, 2), t(2, 2), train_t(3, 2), x(2, 3), tstar(1, 2), test_pair(1, 2);
    VectorXd theta_x = VectorXd::Zero(4);
    VectorXd theta_t = VectorXd::Zero(3);
    VectorXd f = VectorXd::Zero(6);
    t(0, 0) = -0.7258;
    t(0, 1) = -1.9623;
    t(1, 0) = -0.3078;
    t(1, 1) = -0.9332;
    x(0, 0) = 2.4582;
    x(0, 1) = -4.0911;
    x(0, 2) = 1.0004;
    x(1, 0) = 6.1426;
    x(1, 1) = -6.3481;
    x(1, 2) = -4.7591;

    train_t << -0.7258, -1.9623,
    -0.3078, -0.9332,
    0.2501, 1.4158;
    pairs << 0, 1;
    test_pair << 0, 1;
    tstar << 0.2501, 1.4168;
    all_pairs(0) = pairs;
    all_pairs(1) = pairs;


    idx_global_1 << 0, 2;
    idx_global_2 << 1, 3;
    idx_global << 0, 1, 2, 3;
    ind_t << 0, 0, 1, 1;
    ind_x << 0, 1, 0, 1;

    g.Approx_CGppe_Laplace( theta_x, theta_t, sigma,
                           t, x, all_pairs, idx_global, idx_global_1, idx_global_2, ind_t, ind_x, M, N);

    g.Predictive_Utility_Distribution(train_t, tstar, N, idx_global );
    dsp(g.Getmustar(), "mustar");
    dsp(g.Getvarstar(), "varstar");
    return 0;

}


int testpredict_CGppe_laplace_fast()
{
    //generating the data naively
    int M = 3;
    int N = 2;
    double sigma = 0.1;
    CGppe g = CGppe(new CovSEard(), new CovSEard());
    TypePair all_pairs(2);
    VectorXd idx_global_1(2), idx_global_2(2), idx_global(4), ind_t(4), ind_x(4);
    MatrixXd pairs(1, 2), t(2, 2), x(2, 3), tstar(1, 2), test_pair(1, 2);
    VectorXd theta_x = VectorXd::Zero(4);
    VectorXd theta_t = VectorXd::Zero(3);
    VectorXd f = VectorXd::Zero(6);
    t(0, 0) = -0.7258;
    t(0, 1) = -1.9623;
    t(1, 0) = -0.3078;
    t(1, 1) = -0.9332;
    x(0, 0) = 2.4582;
    x(0, 1) = -4.0911;
    x(0, 2) = 1.0004;
    x(1, 0) = 6.1426;
    x(1, 1) = -6.3481;
    x(1, 2) = -4.7591;
    pairs << 0, 1;
    test_pair << 0, 1;
    tstar << 0.2501, 1.4168;
    all_pairs(0) = pairs;
    all_pairs(1) = pairs;


    idx_global_1 << 0, 2;
    idx_global_2 << 1, 3;
    idx_global << 0, 1, 2, 3;
    ind_t << 0, 0, 1, 1;
    ind_x << 0, 1, 0, 1;

    g.Approx_CGppe_Laplace( theta_x, theta_t, sigma,
                           t, x, all_pairs, idx_global, idx_global_1, idx_global_2, ind_t, ind_x, M, N);
             dsp("hello1");
	     dsp(sigma,"sigma");
	    dsp(t,"train_t");
	     dsp(x,"x");
	     dsp(idx_global,"idx_global");
	     dsp(ind_t,"ind_t");
	     dsp(ind_x,"ind_x");
	     dsp(tstar,"test_t");
	     dsp(test_pair,"pair");

    g.Predict_CGppe_Laplace(sigma, t, x, idx_global, ind_t, ind_x, tstar, test_pair);
    return 0;

}







// Some utility routines

void InitMatrix(MatrixXd & mat)
{
    mat(0, 0) = 1.9546;
    mat(0, 1) = 1.5274;
    mat(1, 0) = -0.8292;
    mat(1, 1) = 0.3836;
    mat(2, 0) = 0.9937;
    mat(2, 1) = -1.5854;
    mat(3, 0) = -1.5110;
    mat(3, 1) = -1.3003;
    mat(4, 0) = -2.3473;
    mat(4, 1) = 1.9326;
    mat(5, 0) = 1.2204;
    mat(5, 1) = - 2.3566;
    mat(6, 0) = 0.0001;
    mat(6, 1) = -0.0505;
    mat(7, 0) = - 0.1004;
    mat(7, 1) = - 1.6604;
    mat(8, 0) = 2.0236;
    mat(8, 1) = 2.3934;
    mat(9, 0) = 0.5493;
    mat(9, 1) = 1.0635;
    mat(10, 0) = 0.5883;
    mat(10, 1) = 0.0024;
    mat(11, 0) = 1.7972;
    mat(11, 1) = -0.1446;
}

//Test approc_CGppe_laplace_fast functions of the CGppe class need some fixing
int testapproc_CGppe_laplace_fast()
{
    clock_t start, end;
    double elapsed;
    start = clock();
    //generating the data naively
    int M = 3;
    int N = 2;
    double sigma = 0.1;
    CGppe g = CGppe(new CovSEard(), new CovSEard());
    TypePair all_pairs(2);
    VectorXd idx_global_1(2), idx_global_2(2), idx_global(4), ind_t(4), ind_x(4);
    MatrixXd pairs(1, 2), t(2, 2), x(2, 3), tstar(1, 2);
    VectorXd theta_x = VectorXd::Zero(4);
    VectorXd theta_t = VectorXd::Zero(3);
    t(0, 0) = -0.7258;
    t(0, 1) = -1.9623;
    t(1, 0) = -0.3078;
    t(1, 1) = -0.9332;
    x(0, 0) = 2.4582;
    x(0, 1) = -4.0911;
    x(0, 2) = 1.0004;
    x(1, 0) = 6.1426;
    x(1, 1) = -6.3481;
    x(1, 2) = -4.7591;
    pairs << 0, 1;
    tstar << 0.2501, 1.4168;
    all_pairs(0) = pairs;
    all_pairs(1) = pairs;


    idx_global_1 << 0, 2;
    idx_global_2 << 1, 3;
    idx_global << 0, 1, 2, 3;
    ind_t << 0, 0, 1, 1;
    ind_x << 0, 1, 0, 1;
    g.Approx_CGppe_Laplace( theta_x, theta_t, sigma,
                           t, x, all_pairs, idx_global, idx_global_1, idx_global_2, ind_t, ind_x, M, N);
    dsp(g.GetW(), "W");
    dsp(g.GetL(), "L");
    dsp(g.GetKinv(), "Kinv");
    dsp(g.Getf(), "f");


    end = clock();
    elapsed = ((double)end - start) / CLOCKS_PER_SEC;
    cout << "Elapsed Time :" << elapsed << endl;
    return 0;

}


//Test in order to choose an adaptate support for all_pairs
int testVectors()
{
    MatrixXd v1(1, 3), v2(3, 3), v3;
    v1 << 1, 2, 3;
    v2 << 9, 8, 7, 6, 5, 4, 3, 2, 1;
    v3 = Kron(v1, v2);
    dsp(v3, "v3");
    return 0;

}


// Test voidfunctions
int testvoidfunctions()
{
    MatrixXd z(4, 4), b(5, 5);
    z << 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 14;
    fliplr(z);
    dsp(z, "z");
    return 0;
}



// Test CovSEard
int testCovSEard()
{
	//declaring the data
	int N,M;
	int Maxiter=10;
	TypePair train_pairs, Oracle;
	double sigma=0.1;
 	VectorXd idx_global_1, idx_global_2, idx_global, ind_t, ind_x;
 	MatrixXd pairs(1, 2), test_t, t, x,idx_pairs, F;
 	VectorXd theta_x;
 	VectorXd theta_t;
 	VectorXd ftrue, ytrue;
 	MatrixXd test_pairs;
 	MatrixXd train_t;

 	//assigning the Data
 	Generate(idx_pairs, t, x, Oracle, train_pairs, F, new CovSEard(), new CovSEard(), theta_t, theta_x, M, N, ftrue,ytrue, test_pairs, test_t 
 			,train_t);
 	//Computing the indexes
 	compute_global_index(idx_global_1, idx_global_2, train_pairs, N);
    unique(idx_global, idx_global_1, idx_global_2);
    ind2sub(ind_x, ind_t, N, M, idx_global);
     int test_user_idx=M-1;
     theta_t.fill(0.49782);
    CovSEard a = CovSEard(theta_t);


    double k;
    k = 5;

//cout<<mafunc.Evaluate(t1,t2)<<endl;
    cout << a.ComputeGrandMatrix(train_t) << endl;
//cout<<mafunc.ComputeGrandMatrix(mat)<<endl;
    return 0;
}


//Test CovLINard
int testCovLINard()
{
    VectorXd t1(2), t2(2), t3(3), t4(2), t5(1);
    MatrixXd mat(12, 2);
    for (int z = 0;z < 2;z++)
    {
        t4(z) = 1 ;
    }
    t5(0) = 12;

    InitMatrix(mat);

    t1(0) = 1;
    t1(1) = 15;
    t2(0) = 1;
    t2(1) = 15;
    CovLINard a = CovLINard(t4);
    //CovNoise b=CovNoise(t5);
//CovSum mafunc=CovSum(new CovSEard,new CovSEard,t4);
    //CovSEard mafunc=CovSEard();

    double k;
    k = 5;
//


//cout<<mafunc.Evaluate(t1,t2)<<endl;
    cout << a.ComputeGrandMatrix(mat) << endl;
    cout << a.GetTheta() << endl;

    //cout<<mafunc.ComputeGrandMatrix(mat)<<endl;
    return 0;
}



//Test CoSEiso

int testCoSEiso()
{
    VectorXd t1(2);
    t1 << 1, 15;

    VectorXd t2(2);
    t2 << 1, 15;

    VectorXd t4(2);
    t4 << 1, 1;

    MatrixXd mat(12, 2);
    InitMatrix(mat);

    CovSEiso a = CovSEiso(t4);

//cout<<a.Evaluate(t1,t2)<<endl;
    cout << a.ComputeGrandMatrix(mat) << endl;
    cout << a.GetTheta() << endl;
    return 0;
}

//Test CovNoise

int testCovNoise()
{
    VectorXd t1(2);
    t1 << 1, 15;

    VectorXd t2(2);
    t2 << 1, 15;

    VectorXd t4(1);
    t4 << 1;

    MatrixXd mat(12, 2);
    InitMatrix(mat);

    CovNoise a = CovNoise(t4);

//cout<<a.Evaluate(t1,t2)<<endl;
    cout << a.ComputeGrandMatrix(mat) << endl;
    cout << a.GetTheta() << endl;
    return 0;
}



//Test CovSum

int testCovSum()
{
    VectorXd t1(2);
    t1 << 1, 15;

    VectorXd t2(2);
    t2 << 1, 15;

    VectorXd t5(1);
    t5(0) = 12;

    VectorXd t3(3);
    t3 << 1, 1, 1;

    VectorXd t4(4);
    t4 << 1, 1, 1, 1;

    MatrixXd mat(12, 2);
    InitMatrix(mat);

    CovSum a = CovSum(new CovSEard, new CovNoise, t4);
    CovSEard b = CovSEard(t3);

    cout << a.ComputeGrandMatrix(mat) << endl << endl << endl << endl;
    cout << b.ComputeGrandMatrix(mat) << endl;

    //cout<<a.GetTheta()<<endl;
    return 0;
}


int main()
{
    //testCovSum();
    //testCovNoise();
    //testCoSEiso();
    //testCovLINard();
    //testCovSEard();
    //testVectors();
    //testvoidfunctions();
    //testpredict_CGppe_laplace_fast();
    //testapproc_CGppe_laplace_fast();
    //testpredictive_utility();
    //testmatrixmultiplication();
    //findvalue();
    //testgendata();
    //testnl();
   //testgradnl();
    //testcovderiv();
   //	testopt();
    //testgradcov();
    //testprediction();
    //testinput();
    //testelicit();
   //testreshape();
   //testaddrows();
   // testgenerate();
    //Prediction();
   //Elicitation();
	//Optimisation_without_derivatives();
	//teststring();
	//Posterior();
	//Optimisation();
	//testmat();
	General();
	

}