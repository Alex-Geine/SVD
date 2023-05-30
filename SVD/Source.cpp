#include "Singular.h"

using namespace std;

//метод якоби

int svd2(int m_m, int n_n, double* a, double* u, double* v, double* sigma)
{
	double thr = 1.E-4f, nul = 1.E-16f;
	int col, rows, i, j, l, k, lort, iter, in, ll, kk;
	double alfa, betta, hamma, eta, t, cos0, sin0, buf, s;
	col = n_n;
	rows = m_m;
	for (i = 0; i < col; i++)
	{
		in = i * col;
		for (j = 0; j < col; j++)
			if (i == j) v[in + j] = 1.;
			else v[in + j] = 0.;
	}
	for (i = 0; i < rows; i++)
	{
		in = i * col;
		for (j = 0; j < col; j++)
		{
			u[in + j] = a[in + j];
		}
	}


	iter = 0;
	while (1)
	{
		lort = 0;
		iter++;
		for (l = 0; l < col - 1; l++)
			for (k = l + 1; k < col; k++)
			{
				alfa = 0.; betta = 0.; hamma = 0.;
				for (i = 0; i < rows; i++)
				{
					in = i * col;
					ll = in + l;
					kk = in + k;
					alfa += u[ll] * u[ll];
					betta += u[kk] * u[kk];
					hamma += u[ll] * u[kk];
				}

				if (sqrt(alfa * betta) < nul)	continue;
				if (fabs(hamma) / sqrt(alfa * betta) < thr) continue;

				lort = 1;
				eta = (betta - alfa) / (2.f * hamma);
				t = (double)((eta / fabs(eta)) / (fabs(eta) + sqrt(1. + eta * eta)));
				cos0 = (double)(1. / sqrt(1. + t * t));
				sin0 = t * cos0;

				for (i = 0; i < rows; i++)
				{
					in = i * col;
					buf = u[in + l] * cos0 - u[in + k] * sin0;
					u[in + k] = u[in + l] * sin0 + u[in + k] * cos0;
					u[in + l] = buf;

					if (i >= col) continue;
					buf = v[in + l] * cos0 - v[in + k] * sin0;
					v[in + k] = v[in + l] * sin0 + v[in + k] * cos0;
					v[in + l] = buf;
				}
			}

		if (!lort) break;
	}

	for (i = 0; i < col; i++)
	{
		s = 0.;
		for (j = 0; j < rows; j++)	s += u[j * col + i] * u[j * col + i];
		s = (double)sqrt(s);
		sigma[i] = s;
		if (s < nul)	continue;
		for (j = 0; j < rows; j++)	u[j * col + i] /= s;
	}
	//======= Sortirovka ==============
	for (i = 0; i < col - 1; i++)
		for (j = i; j < col; j++)
			if (sigma[i] < sigma[j])
			{
				s = sigma[i]; sigma[i] = sigma[j]; sigma[j] = s;
				for (k = 0; k < rows; k++)
				{
					s = u[i + k * col]; u[i + k * col] = u[j + k * col]; u[j + k * col] = s;
				}
				for (k = 0; k < col; k++)
				{
					s = v[i + k * col]; v[i + k * col] = v[j + k * col]; v[j + k * col] = s;
				}
			}

	return iter;
}

//m_m - кол-во строк автокоррел€ционной матрицы
//n_n - кол-во столбцов автокоррел€ционной матрицы
//a - автокоррел€ционна€ матрица
//u - лева€ сингул€рна€ матрица
//v - права€ сингул€рна€ матрица
//sigma - матрица с сингул€рными числами

int Svd(int m_m, int n_n, double* a, double* u, double* v, double* sigma)
{

	float thr = 0.000001f; //точность

	int col = n_n, //столбцы
		rows = m_m, //строки
		i,
		j,
		l,
		k,
		lort,
		iter = 0,
		in, //i * n
		ll,
		kk;

	double alfa,
		betta,
		hamma,
		eta,
		t,
		cos0,
		sin0,
		buf,
		s,
		r;

	//заполн€етс€ матрица V, как единична€
	for (i = 0; i < rows; i++) {
		in = i * rows;

		for (j = 0; j < rows; j++)
			if (i == j) v[in + j] = 1.;
			else v[in + j] = 0.;
	}
	
	cout << 'u' << endl;
	//матрица U = A
	for (i = 0; i < rows; i++) {
		in = i * col;

		for (j = 0; j < col; j++) {
			int id = in + j;
			u[in + j] = a[in + j];
			cout << id << ": " << u[in + j] << endl;
		}
			
	}
	//a
	for (int i = 0; i < rows; i++) {
		cout << "| ";
		for (int j = 0; j < col; j++) {
			cout << u[i * col + j] << " ";
		}
		cout << '|' << endl;
	}
	cout << 'v' << endl;
	while (1) {
		lort = 0;
		iter++;


		for (l = 0; l < rows - 1; l++)
			//перебор по строкам начина€ со второй
			for (k = l + 1; k < rows; k++) {
				alfa = 0.;
				betta = 0.;
				hamma = 0.;

				for (int h = 0; h < rows * col + 10; h++) {
					cout << h << ": " << u[h] << endl;
				}
				//перебор столбцов в строке
				for (i = 0; i < col; i++) {
					in = i * rows; // номер первого элеметна в строке
					ll = in + l;
					kk = in + k;
					alfa += u[ll] * u[ll]; //вычисл€етс€ (A * Aт)ll
					betta += u[kk] * u[kk]; //вычисл€етс€ (A * Aт)kk
					hamma += u[ll] * u[kk]; //вычисл€етс€ (A * Aт)lk
				}

				if (sqrt(alfa * betta) < 1.e-10)	continue;
				if (fabs(hamma) / sqrt(alfa * betta) < thr)	continue;

				lort = 1;
				eta = (betta - alfa) / (2.f * hamma);
				t = (eta / fabs(eta)) / (fabs(eta) + (float)sqrt(1.f + eta * eta));
				cos0 = 1.f / (float)sqrt(1.f + t * t);
				sin0 = t * cos0;

				//перебор столбцов
				for (i = 0; i < col; i++) {
					in = i * rows;
					buf = u[in + l] * cos0 - u[in + k] * sin0;
					u[in + k] = u[in + l] * sin0 + u[in + k] * cos0;
					u[in + l] = buf;

					if (i >= rows) continue;
					buf = v[in + l] * cos0 - v[in + k] * sin0;
					v[in + k] = v[in + l] * sin0 + v[in + k] * cos0;
					v[in + l] = buf;
				}
			}

		if (!lort) break;
	}

	//перебор по строкам
	for (i = 0; i < rows; i++) {
		s = 0.;

		for (j = 0; j < col; j++)
			s += u[j * rows + i] * u[j * rows + i];

		s = (float)sqrt(s);
		sigma[i] = s;
		if (s < 1.e-10)	continue;

		for (j = 0; j < col; j++)
			u[j * rows + i] = u[j * rows + i] / s;
	}

	//сортировка
	for (i = 0; i < rows - 1; i++)
		for (j = i; j < rows; j++)
			if (sigma[i] < sigma[j])
			{
				r = sigma[i];
				sigma[i] = sigma[j];
				sigma[j] = r;

				for (k = 0; k < col; k++) {
					r = u[i + k * rows];
					u[i + k * rows] = u[j + k * rows];
					u[j + k * rows] = r;
				}

				for (k = 0; k < rows; k++) {
					r = v[i + k * rows];
					v[i + k * rows] = v[j + k * rows];
					v[j + k * rows] = r;
				}
			}

	return iter;
}

void Show(double* A, int rows, int col) {
	for (int i = 0; i < rows; i++) {
		cout << "| ";
		for (int j = 0; j < col; j++) {
			cout << A[i * col + j] << " ";
		}
		cout << '|' << endl;
	}
}

void Mult(double* U, double* S, double* V, int rows, int col) {
	//трансонирование матрицы V
	double* Vt = new double[col * col];




}

void main() {
	setlocale(LC_ALL, "rus");
	for (;;) {


		int rows, col;
		bool info;
		cout << "¬ведите размерность матрицы: ";
		cin >> rows			//количество строк
			>> col;			//количество столбцоы
		cout<< endl << "¬ывод всех данных(1 - да, 0 - нет): ";
		cin >> info;			//количество строк
		srand(time(0));
		
		Matrix 
			A(rows, col), //матрица формировани€
			X(col, 1),	//вектор неизвестных
			Xres,	//вектор восстановленных неизвестных
			B = A.Mult(X),	//вектор известных
			U, V, Sigma;	//сингул€рное разложение

		A.SVD(U, V, Sigma);		
		

		if (info) {
			cout << 'A' << endl;
			//A.Show();		

			cout << endl << 'U' << endl;
			U.Show();


			cout << endl << 'V' << endl;
			V.Show();

			cout << endl << 'S' << endl;
			Sigma.Show();

			cout << endl << 'x' << endl;
			X.Show();

			cout << endl << 'b' << endl;
			B.Show();
		}
		
		
		//нахождение вектора х
		Matrix Sigmabuf(Sigma.getRows(), Sigma.getCol(), Sigma.conToM());
		
		//Sigma.delRubb();
		Sigma.Trans();
		//cout << endl << 'S' << endl;
		//Sigma.Show();

		Sigma.reverse();
		U.Trans();
		Xres = V.Mult(Sigma.Mult(U.Mult(B)));
		

		if (info) {
			cout << endl << 'x' << endl;
			X.Show();

			cout << endl << 'xres' << endl;
			Xres.Show();
		}

		cout << endl << "Left Error: " << X.lError(Xres) << endl;
		cout << endl << "Right Error: " << Xres.rError(A, B) << endl;
		cout << endl << "Obuslovlennost: " << Sigmabuf.obusl() << endl;
	}
}

