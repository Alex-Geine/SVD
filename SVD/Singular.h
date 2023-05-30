#pragma once

#include <iostream>
#include <time.h>
#include <algorithm>

using namespace std;

//класс для работы с матрицами
class Matrix {
private:
	int
		rows,
		col;
	double** M; //память под матрицу


public:
	//конструкторы
	Matrix();
	Matrix(int rows, int col);
	Matrix(int rows, int col, double* buf);

	//деструктор
	~Matrix();

	//выводит матрицу в консоль
	void Show();

	//get
	int getRows();
	int getCol();
	double** getM();

	//Multiply
	Matrix Mult(Matrix B);

	//транспонирование
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