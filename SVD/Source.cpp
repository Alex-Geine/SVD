#include "Singular.h"
#include <iostream>

using namespace std;



void main() {
	setlocale(LC_ALL, "rus");
	for (;;) {

		try {


			int rows, col;
			bool info;
			cout << "Введите размерность матрицы: ";
			cin >> rows			//количество строк
				>> col;			//количество столбцоы
			cout << endl << "Вывод всех данных(1 - да, 0 - нет): ";
			cin >> info;			//количество строк
			srand(time(0));

			Matrix
				A(rows, col), //матрица формирования
				X(col, 1),	//вектор неизвестных
				Xres,	//вектор восстановленных неизвестных
				B = A.Mult(X),	//вектор известных
				U, V, Sigma;	//сингулярное разложение

			A.SVD(U, V, Sigma);


			if (info) {
				cout << 'A' << endl;
				A.Show();

				cout << endl << 'U' << endl;
				U.Show();


				cout << endl << 'V' << endl;
				V.Show();

				cout << endl << 'S' << endl;
				Sigma.Show();

				cout << endl << 'b' << endl;
				B.Show();
			}


			//нахождение вектора х
			Matrix Sigmabuf(Sigma.getRows(), Sigma.getCol(), Sigma.conToM());


			Sigma.Trans();
			Sigma.reverse();
			U.Trans();
			Xres = V.Mult(Sigma.Mult(U.Mult(B)));


			if (info) {
				cout << endl << 'x' << endl;
				X.Show();

				cout << endl << "xres" << endl;
				Xres.Show();
			}

			cout << endl << "Left Error: " << X.lError(Xres) << endl;
			cout << endl << "Right Error: " << Xres.rError(A, B) << endl;
			cout << endl << "Obuslovlennost: " << Sigmabuf.obusl() << endl;
		}
		catch(...) {

		}
	}
}

