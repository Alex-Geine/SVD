#pragma once

#include <iostream>
#include <time.h>
#include <algorithm>

using namespace std;

//����� ��� ������ � ���������
class Matrix {
private:
	int
		rows,
		col;
	double** M; //������ ��� �������


public:
	//������������
	Matrix();
	Matrix(int rows, int col);
	Matrix(int rows, int col, double* buf);

	//����������
	~Matrix();

	//������� ������� � �������
	void Show();

	//get
	int getRows();
	int getCol();
	double** getM();

	//Multiply
	Matrix Mult(Matrix B);

	//����������������
	void Trans();

	//svd
	void SVD(Matrix &U, Matrix &V, Matrix &Sigma);

	//Convert
	double* conToM();
	double* conToM(bool flag);
	void conToN(double* buf);
	void conToN(double* buf, bool flag);

	//delete rubbish in U matrix
	void delRubb();

	//reverse Sigma
	void reverse();

	//left error
	double lError(Matrix m);
	//right error
	double rError(Matrix A, Matrix B);
	//obusl
	double obusl();
};