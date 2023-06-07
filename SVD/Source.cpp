#include "Singular.h"
#include <iostream>

using namespace std;



void main() {
	setlocale(LC_ALL, "rus");
	for (;;) {

		try {


			int rows, col;
			bool info;
			cout << "������� ����������� �������: ";
			cin >> rows			//���������� �����
				>> col;			//���������� ��������
			cout << endl << "����� ���� ������(1 - ��, 0 - ���): ";
			cin >> info;			//���������� �����
			srand(time(0));

			Matrix
				A(rows, col), //������� ������������
				X(col, 1),	//������ �����������
				Xres,	//������ ��������������� �����������
				B = A.Mult(X),	//������ ���������
				U, V, Sigma;	//����������� ����������

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


			//���������� ������� �
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

