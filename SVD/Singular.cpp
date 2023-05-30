#include "Singular.h"

//конструкторы
Matrix::Matrix(int rows, int col) {
	this->rows = rows;
	this->col = col;

	M = new double*[rows];
	
	for (int i = 0; i < rows; i++) {
		M[i] = new double[col];
		for (int j = 0; j < col; j++) {
			M[i][j] = 1 + (double)rand() / RAND_MAX;			
		}
	}
}

Matrix::Matrix() {

}


Matrix::Matrix(int rows, int col, double* buf) {
	this->rows = rows;
	this->col = col;
	M = new double* [rows];
	//перегон матрицы из одномерного массива
	for (int i = 0; i < rows; i++) {
		M[i] = new double[col];
		for (int j = 0; j < col; j++) {
			int id = i * col + j;
			M[i][j] = buf[id];
		}
	}		
}

//деструктор
Matrix::~Matrix() {
	
}

//выводит матрицу в консоль
void Matrix::Show() {
	for (int i = 0; i < rows; i++) {
		cout << "| ";
		for (int j = 0; j < col; j++) {
			cout << M[i][j] << " ";
		}
		cout << '|' << endl;
	}
}

//get
int Matrix::getRows() {
	return rows;
}

int Matrix::getCol() {
	return col;
}

double** Matrix::getM() {
	return M;
}

//Multiply
Matrix Matrix::Mult(Matrix B) {
	int rowsB = B.getRows();
	int colB = B.getCol();

	if (col > rowsB) {
		cout << "Ошибка перемножения!" << endl;
		return *this;
	}

	Matrix res;
	res.col = colB;
	res.rows = rows;


	res.M = new double *[rows];
	for (int i = 0; i < rows; i++) {
		res.M[i] = new double[colB];
	}

	//m - row, n - col
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < colB; j++) {
			res.M[i][j] = 0;
			for (int k = 0; k < rowsB; k++) {
				if(col-1 >= k)
					res.M[i][j] += M[i][k] * B.M[k][j];
			}
		}
	}	
	return res;
}

//транспонирование
void Matrix::Trans() {
	double** A = new double*[col];


	for (int i = 0; i < col; i++) {
		A[i] = new double[rows];
		for (int j = 0; j < rows; j++) 
			A[i][j] = M[j][i];
	}
	int buf = col;
	col = rows;
	rows = buf;
	M = A;
}

//метод Якоби

//m_m - кол-во строк автокорреляционной матрицы
//n_n - кол-во столбцов автокорреляционной матрицы
//a - автокорреляционная матрица
//u - левая сингулярная матрица
//v - правая сингулярная матрица
//sigma - матрица с сингулярными числами


void Matrix::SVD(Matrix &U, Matrix &V, Matrix &Sigma) {
		//если матрица А переопределенная
		bool transp = false;
		if (this->rows > this->col) {
			this->Trans();	
			transp = true;
		}

		float thr = 0.000001f; //точность

		int Col = this->col, 
			Rows = this->rows, 
			i,
			j,
			l,
			k,
			lort,
			iter = 0,
			in,
			ll,
			kk;		
		
		float alfa,
			betta,
			hamma,
			eta,
			t,
			cos0,
			sin0,
			buf,
			s,
			r;

		V.rows = Col;
		V.col = Col;
		U.col = Col;
		U.rows = Rows;
		Sigma.col = Col;
		Sigma.rows = Rows;

		double
			* a = this->conToM(),
			* v = V.conToM(0),
			* u = U.conToM(0),
			* sigma = Sigma.conToM(0);

		//заполняется матрица V, как единичная
		for (i = 0; i < Col; i++) {
			in = i * Col;

			for (j = 0; j < Col; j++)
				if (i == j) v[in + j] = 1.;
				else v[in + j] = 0.;
		}

		//матрица U = A
		for (i = 0; i < Rows; i++) {
			in = i * Col;

			for (j = 0; j < Col; j++)
				u[in + j] = a[in + j];
		}

		while (1) {
			lort = 0;
			iter++;


			for (l = 0; l < Col - 1; l++)
				//перебор по строкам начиная со второй
				for (k = l + 1; k < Col; k++) {
					alfa = 0.;
					betta = 0.;
					hamma = 0.;

					//перебор столбцов в строке
					for (i = 0; i < Rows; i++) {
						in = i * Col;
						ll = in + l;
						kk = in + k;
						alfa += u[ll] * u[ll];
						betta += u[kk] * u[kk];
						hamma += u[ll] * u[kk];
					}

					if (sqrt(alfa * betta) < 1.e-10)	continue;
					if (fabs(hamma) / sqrt(alfa * betta) < thr)	continue;

					lort = 1;
					eta = (betta - alfa) / (2.f * hamma);
					t = (eta / fabs(eta)) / (fabs(eta) + (float)sqrt(1.f + eta * eta));
					cos0 = 1.f / (float)sqrt(1.f + t * t);
					sin0 = t * cos0;

					for (i = 0; i < Rows; i++)
					{
						in = i * Col;
						buf = u[in + l] * cos0 - u[in + k] * sin0;
						u[in + k] = u[in + l] * sin0 + u[in + k] * cos0;
						u[in + l] = buf;
					}
					for (i = 0; i < Col; i++)
					{
						in = i * Col;
						if (i >= Col) continue;
						buf = v[in + l] * cos0 - v[in + k] * sin0;
						v[in + k] = v[in + l] * sin0 + v[in + k] * cos0;
						v[in + l] = buf;
					}
				}

			if (!lort) break;
		}

		//перебор по строкам
		for (i = 0; i < Col; i++) {
			s = 0.;

			for (j = 0; j < Rows; j++)
				s += u[j * Col + i] * u[j * Col + i];

			s = (float)sqrt(s);
			sigma[i] = s;
			if (s < 1.e-10)	continue;

			for (j = 0; j < Rows; j++)
				u[j * Col + i] = u[j * Col + i] / s;
		}

		//сортировка
		for (i = 0; i < Col - 1; i++)
			for (j = i; j < Col; j++)
				if (sigma[i] < sigma[j])
				{
					r = sigma[i];
					sigma[i] = sigma[j];
					sigma[j] = r;

					for (k = 0; k < Rows; k++) {
						r = u[i + k * Col];
						u[i + k * Col] = u[j + k * Col];
						u[j + k * Col] = r;
					}

					for (k = 0; k < Col; k++) {
						r = v[i + k * Col];
						v[i + k * Col] = v[j + k * Col];
						v[j + k * Col] = r;
					}
				}		

	V.conToN(v);	
	U.conToN(u);	
	Sigma.conToN(sigma, 1);

	if (Rows < Col)
		U.delRubb();

	if (transp) {
		this->Trans();
		Matrix buf = U;
		U = V;
		V = buf;
		Sigma.Trans();
	}
}

//переводит из двойного массива в одинарный
double* Matrix::conToM() {
	double* res = new double[rows * col];
	for (int i = 0; i < rows; i++)
			for (int j = 0; j < col; j++) {
				int id = i * col + j;
				res[id] = M[i][j];
			}
	return res;
}

double* Matrix::conToM(bool flag) {
	double* res = new double[rows * col];	
	return res;
}

//переводит из одинарного в двойной
void Matrix::conToN(double* buf) {
	
	double** A = new double* [rows];
	for (int i = 0; i < rows; i++) {
		A[i] = new double[col];
		for (int j = 0; j < col; j++) {
			int id = i * col + j;
			A[i][j] = buf[id];
		}
	}
	M = A;
}

//переводит из одинарного в двойной диоганальную матрицу
void Matrix::conToN(double* buf, bool flag) {

	double** A = new double* [rows];
	for (int i = 0; i < rows; i++) {
		A[i] = new double[col];
		for (int j = 0; j < col; j++) {
			if (i == j)
				A[i][j] = buf[i];
			else
				A[i][j] = 0;		
		}
	}
	M = A;
}

//удаляет ненужные данные из матрицы U
void Matrix::delRubb() {
	int N;
	int rows = this->rows;
	int col = this->col;

	if (rows > col) {
		this->rows = col;
		N = col;
	}
	else {
		this->col = rows;
		N = rows;
	}
	
	double** A = new double* [N];
	//заполняется матрица V, как единичная
	for (int i = 0; i < N; i++) {
		A[i] = new double[N];		
		for (int j = 0; j < N; j++)
			A[i][j] = this->M[i][j];		
	}
	this->M = A;


	//double* res = new double(rows * rows);

	////заполняется матрица V, как единичная
	//for (int i = 0; i < N; i++) {
	//	int in = i * rows;
	//	for (int j = 0; j < rows; j++)
	//		res[in + j] = this->M[i][j];
	//}

	//this->conToN(res);	
}

//reverse Sigma
void Matrix::reverse() {
	for (int i = 0; i < rows; i++) {		
		for (int j = 0; j < rows; j++)
			if(i == j)
				this->M[i][j] = 1 / this->M[i][j];
	}
}

//left error
double Matrix::lError(Matrix Xres) {
	double res = 0;
	for (int i = 0; i < this->rows; i++)
		res += abs(this->M[i][0] - Xres.M[i][0]);
	return res;
}
//right error
double Matrix::rError(Matrix A, Matrix B) {
	double res = 0;

	Matrix buf = A.Mult(*this);
	for (int i = 0; i < B.rows; i++)
		res += abs(*B.M[i] - *buf.M[i]);
	return res;
}
//obusl
double Matrix::obusl() {
	int N;
	if (this->col > this->rows)
		N = this->rows;
	else
		N = this->col;

	//double min = this->M[N - 1][N - 1];
	//double max = this->M[0][0];
	//cout << "Min: " << this->M[N - 1][N - 1] << " Max: " << this->M[0][0] << endl;

	//this->Show();
	return this->M[0][0] / this->M[N - 1][N - 1];
}