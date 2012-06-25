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

//Test approc_gppe_laplace_fast functions of the Gppe class need some fixing
int main()
{	
	//generating the data naively
	int M=3;
	int N=2;
	double sigma=0.1;
	Gppe g=Gppe();
	TypePair all_pairs(2);
	VectorXd idx_global_1(2),idx_global_2(2),idx_global(4),ind_t(4),ind_x(4),theta_x(4),theta_t(3),f(6);
	MatrixXd pairs(1,2),t(2,2),x(2,3);
	t(0,0)=-0.7258;t(0,1)=-1.9623;t(1,0)=-0.3078;t(1,1)=-0.9332;
	 x(0,0)=2.4582;x(0,1)=-4.0911;x(0,2)=1.0004;
    x(1,0)=6.1426;x(1,1)=-6.3481;x(1,2)=-4.7591;
	pairs<<1,2;
	all_pairs(0)=pairs;
	all_pairs(1)=pairs;	
	
	
	idx_global_1<<0,2;
	idx_global_2<<1,3;
	idx_global<<0,1,2,3;
	ind_t<<0,0,1,1;
	ind_x<<0,1,0,1;

	g.Approx_Gppe_Laplace(new CovSEard,new CovSEard, theta_x, theta_t, sigma,
	t, x, all_pairs, idx_global, idx_global_1, idx_global_2, ind_t, ind_x, M, N);
	//cout<< g.log_likelihood(f, sigma, all_pairs, idx_global_1, idx_global_2, N);
	//cout<<g.GetKx()<<endl;
	return 0;

}









//Test in order to choose an adaptate support for all_pairs
/*
int main()
{
	TypePair mat(2);
	int n=10;
	VectorXd v1(n),v2(n),v3(3);
	v1<<0,1,2,3,4,5,6,7,8,9;
	v2<<0,1,2,3,4,5,6,7,8,9;
	v3<<2,4,6;
	VectorXd v(2);
	v<<n,n;
	v1=2*v1;
	cout<<v1<<endl;
	return 0;

}
*/










// Test voidfunctions
/*
int main()
{
  VectorXd t1(2),t2(2),t3(3),t4(3),t5(1);
  MatrixXd mat(12,2);
  for(int z=0;z<3;z++)
  	{
  		t4(z)=1	;
  	}
  	t5(0)=12;
     mat(0,0)=1.9546;  mat(0,1)=1.5274;
 mat(1,0)=-0.8292; mat(1,1)=0.3836;
  mat(2,0)=0.9937;mat(2,1)=-1.5854;
 mat(3,0)=-1.5110;  mat(3,1)=-1.3003;
mat(4,0)=-2.3473;  mat(4,1)=1.9326;
  mat(5,0)=1.2204;  mat(5,1)=-2.3566;
   mat(6,0)=0.0001;   mat(6,1)=-0.0505;
 mat(7,0)=-0.1004; mat(7,1)=-1.6604;
  mat(8,0)=2.0236;  mat(8,1)=2.3934;
  mat(9,0)=0.5493; mat(9,1)=1.0635;
 mat(10,0)=0.5883;  mat(10,1)=0.0024;
   mat(11,0)=1.7972; mat(11,1)=-0.1446;

  t1(0)=1;
  t1(1)=15;
  t2(0)=1;
  t2(1)=15;
  CovSEard a=CovSEard();
  //CovNoise b=CovNoise(t5);
 //CovSum mafunc=CovSum(new CovSEard,new CovSEard,t4);
   //CovSEard mafunc=CovSEard();
  int k;
  k=5;
 
// a.add(k);
 //cout<<mafunc.Evaluate(t1,t2)<<endl;
  //cout<<a.ComputeGrandMatrix(mat)<<endl;
 //cout<<mafunc.ComputeGrandMatrix(mat)<<endl;
  return 0;
}*/








// Test CovSEard
/*
int main()
{
  VectorXd t1(2),t2(2),t3(3),t4(3),t5(1);
  MatrixXd mat(12,2);
  for(int z=0;z<3;z++)
  	{
  		t4(z)=1	;
  	}
  	t5(0)=12;
     mat(0,0)=1.9546;  mat(0,1)=1.5274;
 mat(1,0)=-0.8292; mat(1,1)=0.3836;
  mat(2,0)=0.9937;mat(2,1)=-1.5854;
 mat(3,0)=-1.5110;  mat(3,1)=-1.3003;
mat(4,0)=-2.3473;  mat(4,1)=1.9326;
  mat(5,0)=1.2204;  mat(5,1)=-2.3566;
   mat(6,0)=0.0001;   mat(6,1)=-0.0505;
 mat(7,0)=-0.1004; mat(7,1)=-1.6604;
  mat(8,0)=2.0236;  mat(8,1)=2.3934;
  mat(9,0)=0.5493; mat(9,1)=1.0635;
 mat(10,0)=0.5883;  mat(10,1)=0.0024;
   mat(11,0)=1.7972; mat(11,1)=-0.1446;

  t1(0)=1;
  t1(1)=15;
  t2(0)=1;
  t2(1)=15;
  CovSEard a=CovSEard(t4);
  //CovNoise b=CovNoise(t5);
 //CovSum mafunc=CovSum(new CovSEard,new CovSEard,t4);
   //CovSEard mafunc=CovSEard();

  double k;
  k=5;
 
 //cout<<mafunc.Evaluate(t1,t2)<<endl;
  cout<<a.ComputeGrandMatrix(mat)<<endl;
 //cout<<mafunc.ComputeGrandMatrix(mat)<<endl;
  return 0;
}*/







//Test CovLINard
/*
int main()
{
  VectorXd t1(2),t2(2),t3(3),t4(2),t5(1);
  MatrixXd mat(12,2);
  for(int z=0;z<2;z++)
  	{
  		t4(z)=1	;
  	}
  	t5(0)=12;
     mat(0,0)=1.9546;  mat(0,1)=1.5274;
 mat(1,0)=-0.8292; mat(1,1)=0.3836;
  mat(2,0)=0.9937;mat(2,1)=-1.5854;
 mat(3,0)=-1.5110;  mat(3,1)=-1.3003;
mat(4,0)=-2.3473;  mat(4,1)=1.9326;
  mat(5,0)=1.2204;  mat(5,1)=-2.3566;
   mat(6,0)=0.0001;   mat(6,1)=-0.0505;
 mat(7,0)=-0.1004; mat(7,1)=-1.6604;
  mat(8,0)=2.0236;  mat(8,1)=2.3934;
  mat(9,0)=0.5493; mat(9,1)=1.0635;
 mat(10,0)=0.5883;  mat(10,1)=0.0024;
   mat(11,0)=1.7972; mat(11,1)=-0.1446;

  t1(0)=1;
  t1(1)=15;
  t2(0)=1;
  t2(1)=15;
  CovLINard a=CovLINard(t4);
  //CovNoise b=CovNoise(t5);
 //CovSum mafunc=CovSum(new CovSEard,new CovSEard,t4);
   //CovSEard mafunc=CovSEard();

  double k;
  k=5;
 //
 
 
 //cout<<mafunc.Evaluate(t1,t2)<<endl;
  cout<<a.ComputeGrandMatrix(mat)<<endl;
  cout<<a.GetTheta()<<endl;

   //cout<<mafunc.ComputeGrandMatrix(mat)<<endl;
  return 0;
}*/







//Test CoSEiso
/*
int main()
{
  VectorXd t1(2),t2(2),t3(3),t4(2),t5(1);
  MatrixXd mat(12,2);
  for(int z=0;z<2;z++)
  	{
  		t4(z)=1	;
  	}
  	t5(0)=12;
     mat(0,0)=1.9546;  mat(0,1)=1.5274;
 mat(1,0)=-0.8292; mat(1,1)=0.3836;
  mat(2,0)=0.9937;mat(2,1)=-1.5854;
 mat(3,0)=-1.5110;  mat(3,1)=-1.3003;
mat(4,0)=-2.3473;  mat(4,1)=1.9326;
  mat(5,0)=1.2204;  mat(5,1)=-2.3566;
   mat(6,0)=0.0001;   mat(6,1)=-0.0505;
 mat(7,0)=-0.1004; mat(7,1)=-1.6604;
  mat(8,0)=2.0236;  mat(8,1)=2.3934;
  mat(9,0)=0.5493; mat(9,1)=1.0635;
 mat(10,0)=0.5883;  mat(10,1)=0.0024;
   mat(11,0)=1.7972; mat(11,1)=-0.1446;

  t1(0)=1;
  t1(1)=15;
  t2(0)=1;
  t2(1)=15;
  CovSEiso a=CovSEiso(t4);
 
 //cout<<a.Evaluate(t1,t2)<<endl;
  cout<<a.ComputeGrandMatrix(mat)<<endl;
  cout<<a.GetTheta()<<endl;
  return 0;
}*/

//Test CovNoise
/*
int main()
{
  VectorXd t1(2),t2(2),t3(3),t4(1),t5(1);
  MatrixXd mat(12,2);
  for(int z=0;z<1;z++)
  	{
  		t4(z)=1	;
  	}
  	t5(0)=12;
     mat(0,0)=1.9546;  mat(0,1)=1.5274;
 mat(1,0)=-0.8292; mat(1,1)=0.3836;
  mat(2,0)=0.9937;mat(2,1)=-1.5854;
 mat(3,0)=-1.5110;  mat(3,1)=-1.3003;
mat(4,0)=-2.3473;  mat(4,1)=1.9326;
  mat(5,0)=1.2204;  mat(5,1)=-2.3566;
   mat(6,0)=0.0001;   mat(6,1)=-0.0505;
 mat(7,0)=-0.1004; mat(7,1)=-1.6604;
  mat(8,0)=2.0236;  mat(8,1)=2.3934;
  mat(9,0)=0.5493; mat(9,1)=1.0635;
 mat(10,0)=0.5883;  mat(10,1)=0.0024;
   mat(11,0)=1.7972; mat(11,1)=-0.1446;

  t1(0)=1;
  t1(1)=15;
  t2(0)=1;
  t2(1)=15;
  CovNoise a=CovNoise(t4);
 
 //cout<<a.Evaluate(t1,t2)<<endl;
  cout<<a.ComputeGrandMatrix(mat)<<endl;
  cout<<a.GetTheta()<<endl;
  return 0;
}*/




//Test CovSum
/*
int main()
{
  VectorXd t1(2),t2(2),t3(3),t4(4),t5(1);
  MatrixXd mat(12,2);
  for(int z=0;z<4;z++)
  	{
  		t4(z)=1	;
  	}
  for(int z=0;z<3;z++)
  	{
  		t3(z)=1	;
  	}
  	t5(0)=12;
     mat(0,0)=1.9546;  mat(0,1)=1.5274;
 mat(1,0)=-0.8292; mat(1,1)=0.3836;
  mat(2,0)=0.9937;mat(2,1)=-1.5854;
 mat(3,0)=-1.5110;  mat(3,1)=-1.3003;
mat(4,0)=-2.3473;  mat(4,1)=1.9326;
  mat(5,0)=1.2204;  mat(5,1)=-2.3566;
   mat(6,0)=0.0001;   mat(6,1)=-0.0505;
 mat(7,0)=-0.1004; mat(7,1)=-1.6604;
  mat(8,0)=2.0236;  mat(8,1)=2.3934;
  mat(9,0)=0.5493; mat(9,1)=1.0635;
 mat(10,0)=0.5883;  mat(10,1)=0.0024;
   mat(11,0)=1.7972; mat(11,1)=-0.1446;

  t1(0)=1;
  t1(1)=15;
  t2(0)=1;
  t2(1)=15;
  CovSum a=CovSum(new CovSEard,new CovNoise,t4);
  CovSEard b=CovSEard(t3);
 
  cout<<a.ComputeGrandMatrix(mat)<<endl<<endl<<endl<<endl;
	cout<<b.ComputeGrandMatrix(mat)<<endl;

  //cout<<a.GetTheta()<<endl;
  return 0;
}*/




