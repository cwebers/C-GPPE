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
#include "Tool.h"





MatrixXd make_query_toydata(TypePair Oracle,int query_idx, int test_idx)
{
	return Oracle(test_idx).row(query_idx);
}

VectorXd get_dWdf(VectorXd all_diag_idx, VectorXd f,VectorXd ind_t,VectorXd ind_x,double sigma, MatrixXd pairs,int  M, int N)
{
	int n=M*N;
	MatrixXd dWdf= MatrixXd::Zero(n,n);
	VectorXd idx_1,idx_2, idx_select;
	
	//We first find the user and the data-point corresponding to this index
	idx_1=pairs.col(0);
	idx_2=pairs.col(1);
	
	// We simply match those corresponding to the required f
	idx_select=find(idx_1,ind_x);
	if(idx_select.rows()==0)
		idx_select=find(idx_2,ind_x);
		
	

	
}



void fliplr(MatrixXd& a)
{
	double inter;
	for(int i=0;i<a.rows();i++)
	{
		for(int j=0;j<((int)(a.cols()/2));j++)
		{
			inter=a(i,j);
			a(i,j)=a(i,a.cols()-1-j);
			a(i,a.cols()-1-j)=inter;
		}
	}

}

void ind2sub(VectorXd& ind_i, VectorXd& ind_j,int dimrow, int dimcol,VectorXd idx )
{
	ind_i=VectorXd::Zero(idx.rows());
	ind_j=VectorXd::Zero(idx.rows());
	int inter;
	
	for(int i=0;i<idx.rows();i++)
	{
	inter=((int)idx(i))%((int)dimrow);
	ind_i(i)=inter;
	inter=idx(i)/dimrow;
	ind_j(i)=inter;
	}
}

void Add(VectorXd& a, double val)
{
	VectorXd inter;
	inter.resize(a.rows()+1);
	inter<<a,val;
	a=inter;
}

VectorXd find(const VectorXd& a, const VectorXd& b)
{
	VectorXd c;
	for(int i=0;i<a.rows();i++)
	{
		if(a(i)==b(i))
			Add(c,i);
	}
		
	return c;
}

int find(const MatrixXd& a,double  val)
{
	int idx=0,i=0,j=0, count=0;
	bool found=false;
	while(i<a.rows()&&found==false)
	{	
		j=0;
		while(j<a.cols()&&found==false)
		{
			if(a(i,j)==val)
			{
				found=true;
				idx=sub2ind(a.rows(),a.cols(), i, j);
			}
			count++;
			j++;
		}
		i++;
	}
	if(found==false)
		idx=INT_MIN;	
	return idx;
}

void unique(VectorXd& a, const VectorXd& b, const VectorXd& c)
{	
	bool check;		
	VectorXd inter;
	a.resize(1,1);
	a<<b(0);
	for(int i=0;i<b.rows();i++)
	{
		check=true;
		for(int j=0;j<a.rows();j++)
		{
			if(b(i)==a(j))
				check=false;
		}
		if(check)
		{
			inter.resize(a.rows()+1,a.cols());
			inter<<a,b(i);
			a=inter;
		}
	
	}
	
	for(int k=0;k<b.rows();k++)
	{
		check=true;
		for(int l=0;l<a.rows();l++)
		{
			if(c(k)==a(l))
				check=false;
		}
		if(check)
		{
			inter.resize(a.rows()+1,a.cols());
			inter<<a,c(k);
			a=inter;
		}
	
	}
			
					
std::sort(a.col(0).data(), a.col(0).data() + a.rows());
}

void dsp(string s)
{
	cout<<endl<<endl<<s<<endl<<endl;
}

void compute_global_index(VectorXd& idx_global_1,VectorXd& idx_global_2,const TypePair& all_pairs,int N)
{
	int M =all_pairs.rows();
    VectorXd inter;

	for(int j=0;j<M;j++)
	{		

			inter.resize(idx_global_1.rows()+ind2global(all_pairs(j).col(0), j, N).rows(),idx_global_1.cols());
			inter<<idx_global_1,ind2global(all_pairs(j).col(0), j, N);
			idx_global_1=inter;
			
			inter.resize(idx_global_2.rows()+ind2global(all_pairs(j).col(1), j, N).rows(),idx_global_2.cols());
			inter<<idx_global_2,ind2global(all_pairs(j).col(1), j, N);
			idx_global_2=inter;
	}
	
}


VectorXd MyNaNMean(MatrixXd a)
{
	VectorXd res=VectorXd::Zero(a.rows());
	int incr;
	double mean=0;
	for (int i=0;i<a.rows();i++)
	{	
		incr=0;
		mean=0;
		for (int j=0;j<a.cols();j++)
		{
				if (a(i,j)!=a(i,j))
					incr++;
				else
					mean+=a(i,j);
		}
		mean=mean/(a.cols()-incr);
		res(i)=mean;
	}
	return res;
}

MatrixXd SetNaN(MatrixXd a)
{
	for (int i=0;i<a.rows();i++)
	{
		for (int j=0;j<a.cols();j++)
		{
				a(i,j)=pow(-1,0.5);
		}
	}
	return a;
	
}


VectorXd GetDiff(VectorXd a,VectorXd b)
{
	VectorXd diff=VectorXd::Zero(a.rows());
	for (int i=0;i<a.rows();i++)
	{
		if(a(i)!=b(i))
			diff(i)=1;
	}
	return diff;
}


void dsp(MatrixXd a,string s)
{
	cout<<endl<<endl<<s<<"	  Size:"<<a.rows()<<" X "<<a.cols()<<endl<<a<<endl<<endl;

}

void dsp(double  a,string s)
{
	cout<<endl<<endl<<s<<" : "<<a<<endl<<endl;
}


MatrixXd GetMatRow(MatrixXd mat,VectorXd t1)
{
	MatrixXd res=MatrixXd::Zero(t1.rows(),mat.cols());
	for(int i=0;i<t1.rows();i++)
	{
		res.row(i)=mat.row(t1(i));
	}
	return res;
}
MatrixXd Kron(MatrixXd mat1, MatrixXd mat2)
{
	MatrixXd mat3=MatrixXd::Zero(mat1.rows()*mat2.rows(), mat1.cols()*mat2.cols());
	for (int i = 0; i < mat1.rows(); i++) 
	{
  		for (int j = 0; j < mat1.cols(); j++) 
  		{	
    		mat3.block(i*mat2.rows(), j*mat2.cols(), mat2.rows(), mat2.cols()) =  mat1(i,j)*mat2;
  		}
	}
	return mat3;
}



VectorXd ind2global(VectorXd vec,int j,int N)
{
 	return vec.array()+j*N;
}

MatrixXd GetMat(MatrixXd mat,VectorXd t1, VectorXd t2)
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

MatrixXd SetMatGenIdx( MatrixXd mat,VectorXd t1, VectorXd t2)
{
	for(int i=0;i<t1.rows();i++)
	{
		mat(t1(i))=t2(i);
	}
	return mat;
}

VectorXd GetMatGenIdx(MatrixXd mat,VectorXd t1)
{
	VectorXd res(t1.rows());
	for(int i=0;i<t1.rows();i++)
	{
		res(i)=mat(t1(i));
	}
	return res;
}

VectorXd GetVec(VectorXd vec,VectorXd t1)
{	
	VectorXd res=VectorXd::Zero(t1.rows());
	for(int i=0;i<t1.rows();i++)
	{
		res(i)=vec(t1(i));
	}
	return res;
}

double normcdf(double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return 0.5*(1.0 + sign*y);
}

VectorXd normcdf(VectorXd x)
{
	for(int i=0;i<x.rows();i++)
	{
		x(i)=normcdf(x(i));
	}
	return x;
}

double normpdf(double x)
{

    return (1./((sqrt(2.*3.14159265358979323846 ))))*exp(-(pow(x,2))/2.);
}

VectorXd normpdf(VectorXd x)
{
	for(int i=0;i<x.rows();i++)
	{
		x(i)=normpdf(x(i));
	}
	return x;
}


VectorXd Get_Cumulative_Val(VectorXd idx, VectorXd val, int n)
{
	VectorXd count=VectorXd::Zero(n);

	for(int i=0;i<val.rows();i++)
	{
		count(idx(i))= count(idx(i)) + val(i); 
	}

	return count;
}

int sub2ind(int dimrow,int dimcol, int row, int col)
{	
	return dimrow*col+row;
}


VectorXd sub2ind(int dimrow,int dimcol, VectorXd setrow,VectorXd setcol)
{
	VectorXd genidx(setrow.rows());
	for(int i=0;i<setrow.rows();i++)
	{
		genidx(i)=sub2ind(dimrow, dimcol,setrow(i),setcol(i));
	}
	return genidx;
}