#include<iostream>
#include<cmath>
#include <Eigen/Eigen>
#include<Eigen/Core>
#include<Eigen/Dense>
using namespace std;
using namespace Eigen;
#define nR 3
#define nT 3
int main()
{
	int i, j, k, min;
	MatrixXcd H(nR, nT), Q(nR, nT), R(nT, nT);
	typedef Matrix<complex<double>, nR, 1> VectornRcd;
	double temp;
	int S[nT];
	VectornRcd x, tR;
	double len;
	typedef Matrix<complex<double>, nT, 1> VectornTcd;
	VectornTcd tT, y, c;
	complex<double> d;
	cout << "Please input the real part of H:" << endl;
	for (i = 0; i < nR; i++)
	{
		for (j = 0; j < nT; j++)
		{
			cin >> temp;
			H.real()(i, j) = temp;
		}
	}
	cout << "Please input the imaginary part of H:" << endl;
	for (i = 0; i < nR; i++)
	{
		for (j = 0; j < nT; j++)
		{
			cin >> temp;
			H.imag()(i, j) = temp;
		}
	}
	cout << "Please input the real part of x:" << endl;
	for (i = 0; i < nR; i++)
	{
		cin >> temp;
		x.real()(i) = temp;
	}
	cout << "Please input the imaginary part of x:" << endl;
	for (i = 0; i < nR; i++)
	{
		cin >> temp;
		x.imag()(i) = temp;
	}
	R.Zero(nT, nT);//R = 0
	Q = H;
	for (i = 0; i < nT; i++)
	{
		S[i] = i;
	}//S=(1,...,nT)
	for (i = 1; i <= nT; i++)
	{
		len = Q.col(i - 1).norm();
		min = i - 1;
		for (j = i; j < nT; j++)
		{
			if (Q.col(j).norm() <= len)
			{
				len = Q.col(j).norm();
				min = j;
			}
		}//compute min = k_i
		tR = Q.col(i - 1);
		Q.col(i - 1) = Q.col(min);
		Q.col(min) = tR;
		tT = R.col(i - 1);
		R.col(i - 1) = R.col(min);
		R.col(min) = tT;
		temp = S[i - 1];
		S[i - 1] = S[min];
		S[min] = temp;//exchange col. i and k_i in Q, R, S
		R(i - 1, i - 1) = Q.col(i - 1).norm();//r_i,i = ||q_i||
		Q.col(i - 1) /= R(i - 1, i - 1);//q_i = q_i / r_i,i
		for (j = i + 1; j <= nT; j++)
		{
			R(i - 1, j - 1) = Q.col(i - 1).adjoint() * Q.col(j - 1);//r_i,l = q_i^H * q_l
			Q.col(j - 1) -= R(i - 1, j - 1) * Q.col(i - 1);//q_l = q_l - r_i,l*q_i
		}
	}//SQRD Algorithm
	y = Q.adjoint() * x;//y = Q^H * x
	c(nT - 1) = y(nT - 1) / R(nT - 1, nT - 1);
	if (abs(c.real()(nT - 1)) > abs(c.imag()(nT - 1)))
	{
		c.real()(nT - 1) = 1;
		c.imag()(nT - 1) = 0;
	}
	else
	{
		c.real()(nT - 1) = 0;
		c.imag()(nT - 1) = 1;
	}//compute c_nT
	for (k = nT - 1; k >= 1; k--)
	{
		d = 0;
		for (i = k + 1; i <= nT; i++)
		{
			d += R(k - 1, i - 1) * c(i - 1);
		}
		d = y(k - 1) - d;
		c(k - 1) = d / R(k - 1, k - 1);
		if (abs(c.real()(k - 1)) > abs(c.imag()(k - 1)))
		{
			c.real()(k - 1) = 1;
			c.imag()(k - 1) = 0;
		}
		else
		{
			c.real()(k - 1) = 0;
			c.imag()(k - 1) = 1;
		}
	}//Signal Detection
	cout << "c =" << endl;
	for (i = 0; i < nT; i++)
	{
		min = 0;
		for (j = 0; j < nT; j++)
		{
			if (S[j] < S[min])
			{
				min = j;
			}
		}
		cout << c(S[min]) << endl;
		S[min] = 101;
	}
}