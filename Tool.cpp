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



VectorXd concatsort(const VectorXd& a,const VectorXd& b)
{
	VectorXd c(a.rows()+b.rows());
	c<<a,b;
	
	std::sort(c.col(0).data(), c.col(0).data() + c.rows());
	return c;	
}

void dspair(TypePair a, string txt)
{
	cout<<txt<<endl;
	for(int i=0;i<a.rows();i++)
	{
		txt="Preférences de lutilisateur num ";
		dsp(a(i),txt);
	}
}
TypePair InputPair(const string& myfile)
{
	int dimpair=0;
	//creating an output flux
	ostringstream oss;
	string info=myfile+"/informations.txt";
	ifstream file(info.c_str());
    if(file)
    {
		file>>dimpair;
    }
    else
    {
        cout << "Couldn't open the file" << endl;
    }
    TypePair res(dimpair);
    for(int i=0;i<dimpair;i++)
    {
    // writing the current number into the flux
    oss << i;
    // Get back the current number into a string 
    info = myfile+"/pair"+oss.str()+".txt";
    res(i)=GetData(info);
    res(i)=res(i).array()-1;
    oss.str("");
    }
    return res;
}


void SetMatRow(MatrixXd& a, VectorXd& t1,MatrixXd& b)
{
	for(int i=0;i<t1.rows();i++)
	{
	a.row(t1(i))=b.row(i);
	}
}


TypePair Gettrain_pairs(TypePair all_pairs)
{
	TypePair train_pairs(all_pairs.rows());
	MatrixXd pairs;
	VectorXd idx;
	for(int i=0;i<train_pairs.rows()-1;i++)
	{
		pairs=all_pairs(i);
		idx=randperm(pairs.rows());
		train_pairs(i)=GetMatRow(pairs,idx);
	}

	return train_pairs;
}


VectorXd randperm(int n)
{
	VectorXd idx(n);
	Vectbool used(n);
	int num=0;
	used.fill(false);
	for(int i=0;i<n;i++)
	{
		srand(time(NULL));
		num=rand() % n;
		while (used(num))
		{
			num=rand() % n;
		}
		idx(i)=num;
		used(num)=true;
	
	}
	return idx;
	
}



MatrixXd combnk(int n)
{
	MatrixXd idx_pairs;
	MatrixXd ln(1,2), inter;
	int z=1;
	for(int i=0;i<n;i++)
	{
		for(int j=z;j<n;j++)
		{
			ln(0,0)=i;
			ln(0,1)=j;
			idx_pairs=MatAdd(idx_pairs,ln);
		}
		z++;
	}
	return idx_pairs;
}


unsigned long BinCoef(int n, int k)
{
	return fact(n)/(fact(k)*fact(n-k));
}


unsigned long fact(int num)
{
    if (num == 0)
        return 1;
    else
        return num * fact(num - 1);
}

void Generate(MatrixXd &idx_pairs,MatrixXd & t,MatrixXd& x,TypePair & Oracle, TypePair & train_pairs, MatrixXd & F, Covfunc *covx, Covfunc *covt, VectorXd &theta_t, VectorXd& theta_x, int &M, int &N,
VectorXd& ftrue_test, VectorXd &ytrue_test, MatrixXd& test_pairs, MatrixXd & test_t, MatrixXd& train_t )
{
	//loading the data
	double EPSILON=0.1405;
	t=GetData("/Users/christopheroustel/Desktop/C-GPPE/t.txt");
	x=GetData("/Users/christopheroustel/Desktop/C-GPPE/x.txt");
	M=t.rows();
	N=x.rows();
	int n=M*N;
	theta_t=VectorXd::Zero(t.cols()+1);
	theta_x=VectorXd::Zero(x.cols()+1);
	covt->SetTheta(theta_t);
	covx->SetTheta(theta_x);
	MatrixXd Kt=covt->ComputeGrandMatrix(t);
	MatrixXd Kx=covx->ComputeGrandMatrix(x);
	MatrixXd K=Kron(Kt,Kx);
	F=GetData("/Users/christopheroustel/Desktop/C-GPPE/F.txt");
	//other way to get F with random numbers
	//LLT<MatrixXd> Kllt;
	//Kllt.compute(K);
	//K=Kllt.matrixL();
	//VectorXd rdnum(n);
	//rdnum.setRandom();
	//VectorXd f=K*rdnum;
	//F=reshape(f,N,M);
	//idx_pairs=combnk(N);
	idx_pairs=GetData("/Users/christopheroustel/Desktop/C-GPPE/idx_pairs.txt");
	idx_pairs=idx_pairs.array()-1;
	MatrixXd Y=(GetMatRow(F,idx_pairs.col(0))-GetMatRow(F,idx_pairs.col(1)));
	for(int i=0;i<Y.rows();i++)
	{
		for(int j=0;j<Y.cols();j++)
		{
			if(Y(i,j)>EPSILON)
				Y(i,j)=1;
			else
				Y(i,j)=0;
		}
	}
	
	//generating the preferences pairs
	Oracle=TypePair(M);
	MatrixXd tmp_pairs;
	VectorXd idx_0;
	MatrixXd inter;
 	//Myself
	//for(int i=0;i<M;i++)
	//{
	//	tmp_pairs=idx_pairs;
	//	idx_0=find(Y.col(i),0);
	//	inter=GetMatRow(idx_pairs,idx_0);

	//	fliplr(inter);
	//	SetMatRow(tmp_pairs,idx_0,inter);
	//Oracle(i)=tmp_pairs;
	//}
		//with matlab's data
	Oracle=InputPair("/Users/christopheroustel/Desktop/C-GPPE/all_pairs");
	//Generating train_pairs myself
	//train_pairs=Gettrain_pairs(Oracle);
	//with matlab's data
	train_pairs=InputPair("/Users/christopheroustel/Desktop/C-GPPE/train_pairs");

int Mtrain = M-1;
int test_idx = M-1;
test_pairs  = Oracle(test_idx);
train_t = GetMatRow(t,Nfirst(Mtrain));
test_t = t.row(test_idx);
ftrue_test = F.col(test_idx);
ytrue_test = Y.col(test_idx);

}

MatrixXd MatAdd(MatrixXd  mat, MatrixXd ln)
{
	//resizing the matrix without erazing the previous data
	mat.conservativeResize(mat.rows()+1,ln.cols());
	//Adding the new line
	for(int z=0;z<ln.cols();z++)
	{
		mat(mat.rows()-1,z)=ln(0,z);
	}
	return mat;
}



MatrixXd reshape(VectorXd f, int ln, int col)
{
	MatrixXd a(ln,col);
	int z=0;
	for(int j=0;j<col;j++)
	{
		for(int i=0;i<ln;i++)
		{
			a(i,j)=f(z);
			z++;
		}
	}
	return a;
}

MatrixXd GetData(const string& myfile)
{
	int dim_line;
	int dim_col;
	MatrixXd a;
	ifstream file(myfile.c_str());

    if(file)
    {
    	file>>dim_line;
    	file>>dim_col;
		a=MatrixXd::Zero(dim_line,dim_col);//now we can create the Matrix
		for(int i=0;i<dim_line;i++)
		{
			for (int j=0;j<dim_col;j++)
			{
				file>>a(i,j);
			}
		}
    }
    else
    {
        cout << "Couldn't open the file" << endl;
    }
    
    return a;
}


VectorXd find(const VectorXd a, double b)
{
	VectorXd res;
	int z=0;
	for (int i=0;i<a.rows();i++)
	{
		if(a(i)==b)
		{
			Add(res,i);
		}
	}
	return res;
}



column_vector EigentoDlib(VectorXd a)
{
    column_vector b(a.rows());
    for (int i = 0;i < a.rows();i++)
    {
        b(i) = a(i);
    }
    return b;
}

VectorXd DlibtoEigen(column_vector a)
{
    VectorXd b(a.size());
    for (int i = 0;i < b.rows();i++)
    {
        b(i) = a(i);
    }
    return b;
}

void GetTheta(VectorXd& theta_t, VectorXd& theta_x, double& sigma , VectorXd& theta, int dim_t, int dim_x)
{
    double dim = theta.rows();
    sigma = theta(dim - 1);
    theta_t = theta.head(dim_t+1);
    theta_x = theta.segment(dim_t+1,dim_x+1);

}

VectorXd concatTheta(const VectorXd &theta_t, const VectorXd &theta_x, double sigma)
{
    VectorXd theta(theta_t.rows() + theta_x.rows() + 1);
    theta << theta_t, theta_x, sigma;
    return theta;
}



VectorXd concatmat(const MatrixXd& a )
{
    VectorXd concat(a.rows()*a.cols());
    int z = 0;
    for (int j = 0;j < a.cols();j++)
    {
        for (int i = 0;i < a.rows();i++)
        {
            concat(z) = a(i, j);
            z++;
        }
    }
    return concat;
}



VectorXd Nfirst(int N)
{
    VectorXd a(N);
    for (int i = 0;i < N;i++)
    {
        a(i) = i;
    }
    return a;
}

VectorXd get_dlogp_dsigma(MatrixXd& dWdsigma, double& dloglike_dsigma, const VectorXd& f, double sigma, const TypePair& all_pairs, int M, int N)
{
    M = all_pairs.rows();
    MatrixXd Dfdsigma = MatrixXd::Zero(N, M);
    VectorXd idx_global_1, idx_global_2, pdf_val, cdf_val, ratio;
    VectorXd idx_1, idx_2, z, val, coef;
    MatrixXd pairs;

    for (int j = 0;j < M;j++)
    {
        if (all_pairs(j).rows() == 0)
            continue;
        pairs = all_pairs(j);
        idx_1 = pairs.col(0);
        idx_2 = pairs.col(1);

        idx_global_1 = ind2global(idx_1, j, N);
        idx_global_2 = ind2global(idx_2, j, N);

        z = (GetVec(f, idx_global_1) - GetVec(f, idx_global_2)) / sigma;
        pdf_val = normpdf(z);
        cdf_val = normcdf(z);

        ratio = pdf_val.array() / cdf_val.array();
        VectorXd inter = (z.array() * (z + ratio).array()).array() - 1;

        val = (1 / sigma) * ratio.array() * inter.array();

        coef = Get_Cumulative_Val(idx_1, val, N);
        coef = coef - Get_Cumulative_Val(idx_2, val, N);
        Dfdsigma.col(j) = coef;
    }
    return concatmat(Dfdsigma);
}


void get_dsigma(MatrixXd& dWdsigma, double& dloglike_dsigma, const VectorXd& f, double sigma, const TypePair& all_pairs, int M, int N)
{
    M = all_pairs.rows();
    VectorXd idx_1, idx_2, idx_global_1, idx_global_2, z, pdf_val, cdf_val;
    VectorXd ratio1, val, ind, ind_trans;
    int n = M * N;
    MatrixXd pairs;
    dWdsigma = MatrixXd::Zero(n, n);

    VectorXd all_diag_idx = sub2ind(n, n, Nfirst(n), Nfirst(n));
    dloglike_dsigma = 0;

    for (int j = 0;j < M;j++)
    {
        if (all_pairs(j).rows() == 0)
            continue;
        pairs = all_pairs(j);
        idx_1 = pairs.col(0);
        idx_2 = pairs.col(1);

        idx_global_1 = ind2global(idx_1, j, N);
        idx_global_2 = ind2global(idx_2, j, N);

        z = (GetVec(f, idx_global_1) - GetVec(f, idx_global_2)) / sigma;
        pdf_val = normpdf(z);
        cdf_val = normcdf(z);

        ratio1 = pdf_val.array() / cdf_val.array();
        VectorXd inter1, inter2, inter3;
        inter2 = (1 - ((z + ratio1).array() * (z + 2 * ratio1).array()).array());
        inter3 = 2 * ratio1.array() * (z + ratio1).array();
        inter1 = z.array() * ratio1.array();

        val = (-1 / pow(sigma, 2)) *( inter1.array() * inter2.array() + inter3.array());
        ind = sub2ind(n, n, idx_global_1, idx_global_2);
        dWdsigma = SetMatGenIdx(dWdsigma, ind, -1 * val);

        ind_trans = sub2ind(n, n, idx_global_2, idx_global_1);
        dWdsigma = SetMatGenIdx(dWdsigma, ind_trans, -1 * val);

        //Now computing the diagonal
        dWdsigma = SetMatGenIdx(dWdsigma, all_diag_idx, GetMatGenIdx(dWdsigma, all_diag_idx) + Get_Cumulative_Val(idx_global_1, val, n));
        dWdsigma = SetMatGenIdx(dWdsigma, all_diag_idx, GetMatGenIdx(dWdsigma, all_diag_idx) + Get_Cumulative_Val(idx_global_2, val, n));

        dloglike_dsigma = dloglike_dsigma - (z.array() * ratio1.array()).sum();

    }

}





void loss_query_toydata(double& loss, const MatrixXd& F, bool& stop, int test_user_idx, int best_item_idx)
{
    VectorXd ftest = F.col(test_user_idx);
    double best_val = ftest.maxCoeff();
    double pred_val = ftest(best_item_idx);
    loss = best_val - pred_val;
    if (pred_val == best_val)
        stop = true;
    else
        stop = false;
}



MatrixXd make_query_toydata(TypePair Oracle, int query_idx, int test_idx)
{
    return Oracle(test_idx).row(query_idx);
}

MatrixXd get_dWdf(VectorXd all_diag_idx, VectorXd f, int ind_t, int ind_x, double sigma, MatrixXd pairs, int  M, int N)
{
    int n = M * N;
    MatrixXd dWdf = MatrixXd::Zero(n, n);
    VectorXd idx_1, idx_2, idx_select,idx_select_1, idx_select_2, idx_global_1, idx_global_2, ind, z, pdf_val, cdf_val, val, coeff, ratio1, ind_trans;

    //We first find the user and the data-point corresponding to this index
    idx_1 = pairs.col(0);
    idx_2 = pairs.col(1);

    // We simply match those corresponding to the required f
    idx_select_1 = find(idx_1, ind_x);
	idx_select_2 = find(idx_2, ind_x);
	idx_select=concatsort(idx_select_1,idx_select_2);

    if (idx_select.rows() == 0)
        return dWdf;

    idx_1 = GetVec(idx_1, idx_select);
    idx_2 = GetVec(idx_2, idx_select);

    coeff=VectorXd::Zero(idx_1.rows());
    coeff.fill(1);
    VectorXd simili = find(idx_2, ind_x);

    //it is negative if f_{o} is on the wrong side of the relationship
    for (int i = 0;i < simili.rows();i++)
    {
        coeff(simili(i)) = -1;
    }

    idx_global_1 = ind2global(idx_1, ind_t, N);
    idx_global_2 = ind2global(idx_2, ind_t, N);
    z = (GetVec(f, idx_global_1) - GetVec(f, idx_global_2)) / sigma;
    pdf_val = normpdf(z);
    cdf_val = normcdf(z);
    ratio1 = pdf_val.array() / cdf_val.array();

    val = (1 / pow(sigma, 3)) * ratio1.array() * (1 - ((z + ratio1).array() * (z + 2 * ratio1).array()).array()).array();

    val = val.array() * coeff.array();

    ind = ind2global(idx_global_1, idx_global_2, n);
    dWdf = SetMatGenIdx(dWdf, ind, -1 * val);

    ind_trans = ind2global(idx_global_2, idx_global_1, n);
    dWdf = SetMatGenIdx(dWdf, ind_trans, -1 * val);


//now taking care of the diagonal

    dWdf = SetMatGenIdx(dWdf, all_diag_idx, GetMatGenIdx(dWdf, all_diag_idx) + Get_Cumulative_Val(idx_global_1, val, n));
    dWdf = SetMatGenIdx(dWdf, all_diag_idx, GetMatGenIdx(dWdf, all_diag_idx) + Get_Cumulative_Val(idx_global_2, val, n));
    return dWdf;

}



void fliplr(MatrixXd& a)
{
    double inter;
    for (int i = 0;i < a.rows();i++)
    {
        for (int j = 0;j < ((int)(a.cols() / 2));j++)
        {
            inter = a(i, j);
            a(i, j) = a(i, a.cols() - 1 - j);
            a(i, a.cols() - 1 - j) = inter;
        }
    }

}

void ind2sub(VectorXd& ind_i, VectorXd& ind_j, int dimrow, int dimcol, VectorXd idx )
{
    ind_i = VectorXd::Zero(idx.rows());
    ind_j = VectorXd::Zero(idx.rows());
    int inter;

    for (int i = 0;i < idx.rows();i++)
    {
        inter = ((int)idx(i)) % ((int)dimrow);
        ind_i(i) = inter;
        inter = idx(i) / dimrow;
        ind_j(i) = inter;
    }
}

void Add(VectorXd& a, double val)
{
    VectorXd inter;
    inter.resize(a.rows() + 1);
    inter << a, val;
    a = inter;
}


void unique(VectorXd& a, const VectorXd& b, const VectorXd& c)
{
    bool check;
    VectorXd inter;
    a.resize(1, 1);
    a << b(0);
    for (int i = 0;i < b.rows();i++)
    {
        check = true;
        for (int j = 0;j < a.rows();j++)
        {
            if (b(i) == a(j))
                check = false;
        }
        if (check)
        {
            inter.resize(a.rows() + 1, a.cols());
            inter << a, b(i);
            a = inter;
        }

    }

    for (int k = 0;k < b.rows();k++)
    {
        check = true;
        for (int l = 0;l < a.rows();l++)
        {
            if (c(k) == a(l))
                check = false;
        }
        if (check)
        {
            inter.resize(a.rows() + 1, a.cols());
            inter << a, c(k);
            a = inter;
        }

    }


    std::sort(a.col(0).data(), a.col(0).data() + a.rows());
}

void dsp(string s)
{
    cout << endl << endl << s << endl << endl;
}

void compute_global_index(VectorXd& idx_global_1, VectorXd& idx_global_2, const TypePair& all_pairs, int N)
{
    int M = all_pairs.rows();
    VectorXd inter;
	idx_global_1.resize(0,0);
	idx_global_2.resize(0,0);
    for (int j = 0;j < M;j++)
    {

        inter.resize(idx_global_1.rows() + ind2global(all_pairs(j).col(0), j, N).rows(), idx_global_1.cols());
        inter << idx_global_1, ind2global(all_pairs(j).col(0), j, N);
        idx_global_1 = inter;

        inter.resize(idx_global_2.rows() + ind2global(all_pairs(j).col(1), j, N).rows(), idx_global_2.cols());
        inter << idx_global_2, ind2global(all_pairs(j).col(1), j, N);
        idx_global_2 = inter;
    }

}


VectorXd GetDiff(VectorXd a, VectorXd b)
{
    VectorXd diff = VectorXd::Zero(a.rows());
    for (int i = 0;i < a.rows();i++)
    {
        if (a(i) != b(i))
            diff(i) = 1;
    }
    return diff;
}


void dsp(MatrixXd a, string s)
{
    cout << s << "   Size:" << a.rows() << " X " << a.cols() << endl << a << endl << endl;

}

void dsp(double  a, string s)
{
    cout << s << " : " << a << endl << endl;
}


MatrixXd GetMatRow(MatrixXd mat, VectorXd t1)
{
    MatrixXd res = MatrixXd::Zero(t1.rows(), mat.cols());
    for (int i = 0;i < t1.rows();i++)
    {
        res.row(i) = mat.row(t1(i));
    }
    return res;
}
MatrixXd Kron(MatrixXd mat1, MatrixXd mat2)
{
    MatrixXd mat3 = MatrixXd::Zero(mat1.rows() * mat2.rows(), mat1.cols() * mat2.cols());
    for (int i = 0; i < mat1.rows(); i++)
    {
        for (int j = 0; j < mat1.cols(); j++)
        {
            mat3.block(i*mat2.rows(), j*mat2.cols(), mat2.rows(), mat2.cols()) =  mat1(i, j) * mat2;
        }
    }
    return mat3;
}



VectorXd ind2global(VectorXd vec, int j, int N)
{
    return vec.array() + j*N;
}

VectorXd ind2global(VectorXd a, VectorXd b, int N)
{
    return a + b*N;
}

MatrixXd GetMat(MatrixXd mat, VectorXd t1, VectorXd t2)
{
    MatrixXd res(t1.rows(), t2.rows());
    for (int i = 0;i < t1.rows();i++)
    {
        for (int j = 0;j < t2.rows();j++)
        {
            res(i, j) = mat(t1(i), t2(j));
        }
    }
    return res;
}

MatrixXd SetMatGenIdx( MatrixXd mat, VectorXd t1, VectorXd t2)
{
    for (int i = 0;i < t1.rows();i++)
    {
        mat(t1(i)) = t2(i);
    }
    return mat;
}

VectorXd GetMatGenIdx(MatrixXd mat, VectorXd t1)
{
    VectorXd res(t1.rows());
    for (int i = 0;i < t1.rows();i++)
    {
        res(i) = mat(t1(i));
    }
    return res;
}

VectorXd GetVec(VectorXd vec, VectorXd t1)
{
    VectorXd res = VectorXd::Zero(t1.rows());
    for (int i = 0;i < t1.rows();i++)
    {
        res(i) = vec(t1(i));
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
    x = fabs(x) / sqrt(2.0);

    // A&S formula 7.1.26
    double t = 1.0 / (1.0 + p * x);
    double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-x * x);

    return 0.5*(1.0 + sign*y);
}

VectorXd normcdf(VectorXd x)
{
    for (int i = 0;i < x.rows();i++)
    {
        x(i) = normcdf(x(i));
    }
    return x;
}

double normpdf(double x)
{

    return (1. / ((sqrt(2.*3.14159265358979323846 ))))*exp(-(pow(x, 2)) / 2.);
}

VectorXd normpdf(VectorXd x)
{
    for (int i = 0;i < x.rows();i++)
    {
        x(i) = normpdf(x(i));
    }
    return x;
}


VectorXd Get_Cumulative_Val(VectorXd idx, VectorXd val, int n)
{
    VectorXd count = VectorXd::Zero(n);

    for (int i = 0;i < val.rows();i++)
    {
        count(idx(i)) = count(idx(i)) + val(i);
    }

    return count;
}

int sub2ind(int dimrow, int dimcol, int row, int col)
{
    return dimrow*col + row;
}


VectorXd sub2ind(int dimrow, int dimcol, VectorXd setrow, VectorXd setcol)
{
    VectorXd genidx(setrow.rows());
    for (int i = 0;i < setrow.rows();i++)
    {
        genidx(i) = sub2ind(dimrow, dimcol, setrow(i), setcol(i));
    }
    return genidx;
}