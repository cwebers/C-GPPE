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
#include "Gppe.h"

int bigapproc_gppe_laplace_fast()
{
	//clock_t start, end;
	//double elapsed;
	//start = clock();
    //generating the data naively
    int M = 4;
    int N = 4;
    double sigma = 0.1;
    Gppe g = Gppe(new CovSEard(), new CovSEard());
    TypePair all_pairs(4);
    VectorXd idx_global_1(18), idx_global_2(18), idx_global(12), ind_t(12), ind_x(12);
    MatrixXd pairs(6, 2),pairs2(6, 2),pairs3(6, 2),pairs4(6, 2), t(3, 2), x(4, 3), tstar(1,2);
    VectorXd theta_x = VectorXd::Zero(4);
    VectorXd theta_t = VectorXd::Zero(3);
    VectorXd f ;
   
   
   t<<   -0.7258 ,  -0.9332,
   		-0.3078 ,   1.4168,
   		 0.2501 ,   0.8194;
   		 
   		 
   	x<<   -4.0911,  0.4715,   -7.3005,
   		-6.3481 ,   4.8933,    5.8177,
    	1.0004 ,   4.3248,    2.3851,
   		-4.7591,   -3.7461,   -6.3772;
    tstar<< 0.2501, 1.4168;
     pairs<< 3,4,4,2,3,2,1,4,1,3,1,2;
     all_pairs(0)=pairs;
    pairs2 << 3,4,4,2,3,2,4,1,3,1,1,2;
    all_pairs(1)=pairs2;
    
    pairs3<<  4,3,4,2,3,2,4,1,3,1,2,1; 
    all_pairs(2)=pairs3;
    pairs4<<  3,4,2,4,2,3,1,4,3,1,2,1;
    all_pairs(3)=pairs4;


    idx_global_1 <<2,0,0,3,2,0,6,7,7,6,4,6,11,10,9,10,11,11;
    idx_global_2 <<3,1,2,1,1,3,4,4,5,7,5,5,10,8,8,9,8,9;
    idx_global << 0, 1, 2, 3,4,5,6,7,8,9,10,11;
    ind_t <<0,0,0,0,1,1,1,1,2,2,2,2;
    ind_x <<0,1,2,3,0,1,2,3,0,1,2,3;
	cout<<"Hello"<<endl;
    g.Approx_Gppe_Laplace( theta_x, theta_t, sigma,
                          t, x, all_pairs, idx_global, idx_global_1, idx_global_2, ind_t, ind_x, M, N);
    
    dsp(g.GetW(),"W");
    dsp(g.GetL(),"L");
    dsp(g.GetKinv(),"Kinv");
    dsp(g.Getf(),"f");
  
  
      // end = clock();
    // elapsed = ((double)end - start) / CLOCKS_PER_SEC;
  //   cout << elapsed << endl;
  return 0;
}

int testpredict_gppe_laplace_fast()
{
    //generating the data naively
    int M = 3;
    int N = 2;
    double sigma = 0.1;
    Gppe g = Gppe(new CovSEard(), new CovSEard());
    TypePair all_pairs(2);
    VectorXd idx_global_1(2), idx_global_2(2), idx_global(4), ind_t(4), ind_x(4);
    MatrixXd pairs(1, 2), t(2, 2), x(2, 3), tstar(1,2),test_pair(1,2);
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
    pairs << 1, 2;
    test_pair << 0, 1;
    tstar<< 0.2501, 1.4168;
    all_pairs(0) = pairs;
    all_pairs(1) = pairs;


    idx_global_1 << 0, 2;
    idx_global_2 << 1, 3;
    idx_global << 0, 1, 2, 3;
    ind_t << 0, 0, 1, 1;
    ind_x << 0, 1, 0, 1;

    g.Approx_Gppe_Laplace( theta_x, theta_t, sigma,
                          t, x, all_pairs, idx_global, idx_global_1, idx_global_2, ind_t, ind_x, M, N);
		
	g.Predict_Gppe_Laplace(sigma, t, x,idx_global, ind_t, ind_x,tstar, test_pair);
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

//Test approc_gppe_laplace_fast functions of the Gppe class need some fixing
int testapproc_gppe_laplace_fast()
{
	//clock_t start, end;
	//double elapsed;
	//start = clock();
    //generating the data naively
    int M = 3;
    int N = 2;
    double sigma = 0.1;
    Gppe g = Gppe(new CovSEard(), new CovSEard());
    TypePair all_pairs(2);
    VectorXd idx_global_1(2), idx_global_2(2), idx_global(4), ind_t(4), ind_x(4);
    MatrixXd pairs(1, 2), t(2, 2), x(2, 3), tstar(1,2);
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
    pairs << 1, 2;
    tstar<< 0.2501, 1.4168;
    all_pairs(0) = pairs;
    all_pairs(1) = pairs;


    idx_global_1 << 0, 2;
    idx_global_2 << 1, 3;
    idx_global << 0, 1, 2, 3;
    ind_t << 0, 0, 1, 1;
    ind_x << 0, 1, 0, 1;
    
    g.Approx_Gppe_Laplace( theta_x, theta_t, sigma,
                          t, x, all_pairs, idx_global, idx_global_1, idx_global_2, ind_t, ind_x, M, N);
    dsp(g.GetW(),"W");
    dsp(g.GetL(),"L");
    dsp(g.GetKinv(),"Kinv");
    dsp(g.Getf(),"f");
  
  
      // end = clock();
    // elapsed = ((double)end - start) / CLOCKS_PER_SEC;
  //   cout << elapsed << endl;
  return 0;

}


//Test in order to choose an adaptate support for all_pairs
int testVectors()
{
    MatrixXd v1(1,3), v2(3,3), v3;
    v1 <<1,2,3;
    v2 <<9,8,7,6,5,4,3,2,1;
    v3=Kron(v1,v2);
     dsp(v3,"v3");
    return 0;

}


// Test voidfunctions
int testvoidfunctions()
{
    VectorXd t1(2), t2(2), t3(3), t4(3), t5(1);
    MatrixXd mat(12, 2);
    for (int z = 0;z < 3;z++)
    {
        t4(z) = 1 ;
    }
    t5(0) = 12;/*
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
mat(11, 1) = -0.1446;*/

    t1(0) = 1;
    t1(1) = 15;
    t2(0) = 1;
    t2(1) = 15;
    // CovSEard a = CovSEard();
    //CovNoise b=CovNoise(t5);
//CovSum mafunc=CovSum(new CovSEard,new CovSEard,t4);
    //CovSEard mafunc=CovSEard();
    int k;
    k = 5;

// a.add(k);
//cout<<mafunc.Evaluate(t1,t2)<<endl;
    //cout<<a.ComputeGrandMatrix(mat)<<endl;
//cout<<mafunc.ComputeGrandMatrix(mat)<<endl;
    return 0;
}



// Test CovSEard
int testCovSEard()
{
    VectorXd t1(2), t2(2), t3(3), t4(3), t5(1);
    MatrixXd mat(12, 2);
    for (int z = 0;z < 3;z++)
    {
        t4(z) = 1 ;
    }
    t5(0) = 12;

    InitMatrix(mat);

    t1(0) = 1;
    t1(1) = 15;
    t2(0) = 1;
    t2(1) = 15;
    CovSEard a = CovSEard(t4);
    //CovNoise b=CovNoise(t5);
//CovSum mafunc=CovSum(new CovSEard,new CovSEard,t4);
    //CovSEard mafunc=CovSEard();

    double k;
    k = 5;

//cout<<mafunc.Evaluate(t1,t2)<<endl;
    cout << a.ComputeGrandMatrix(mat) << endl;
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
    //testpredict_gppe_laplace_fast();
    //testapproc_gppe_laplace_fast();
	bigapproc_gppe_laplace_fast();

}

