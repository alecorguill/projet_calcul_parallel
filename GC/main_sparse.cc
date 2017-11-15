#include <iostream>
#include<complex>
#include "Dense"
#include "Sparse"
using namespace Eigen;
using namespace std;

int main()
{
	SparseMatrix<double> Id(5,5), O(15,15), S(10,10), T(10,10);   // (30,30) sparse matrix
	SparseVector<double> x(10);
	Matrix<double, 10, 10> A;

	cout << "-------------------------------" << endl;
	O.setZero();    // Matrix with zeros
	cout << "O (setZero) = " << endl;
	cout << O << endl;
	cout << "-------------------------------" << endl;
	Id.setIdentity();    // Matrix with zeros
	cout << "Id (setIdentity) = " << endl;
	cout << Id << endl;
	cout << "-------------------------------" << endl;

	vector<Triplet<double>> triplets;
	for (int i=0; i<S.rows(); ++i) {
		triplets.push_back({i,i,1.});
		if (i > 0)
		triplets.push_back({i,i-1,-2.});
		if (i < S.rows()-1)
		triplets.push_back({i,i+1,2.});
	}
	S.setFromTriplets(triplets.begin(), triplets.end());

	triplets.clear();
	for (int i=0; i<T.rows(); ++i) {
		triplets.push_back({i,i,2.});
		if (i > 1)
			triplets.push_back({i,i-2,-4.});
		if (i < T.rows()-2)
			triplets.push_back({i,i+2,4.});
		if (i > 3)
			triplets.push_back({i,i-4,-8.});
		if (i < T.rows()-4)
			triplets.push_back({i,i+4,8.});
	}
	T.setFromTriplets(triplets.begin(), triplets.end());

	vector<int> index_x = {0, 2, 9};
	for (int i = 0 ; i < index_x.size() ; ++i)
		x.coeffRef(index_x[i]) = 2.*i;

	A = MatrixXd::Random(A.rows(),A.cols());

	cout << "-------------------------------" << endl;
	cout << "S (build with triplets) = " << endl;
	cout << S << endl;
	cout << "-------------------------------" << endl;
	cout << "T (build with triplets) = " << endl;
	cout << T << endl;
	cout << "-------------------------------" << endl;
	cout << "x (build with coeffRef) = " << endl;
	cout << x << endl;
	cout << "-------------------------------" << endl;
	cout << "A (random DENSE matrix) = " << endl;
	cout << A << endl;
	cout << "-------------------------------" << endl;

	cout << "-------------------------------" << endl;
	cout << "Adjoint, transpose " << endl;
	cout << "-------------------------------" << endl;
	cout << "S^T = " << endl;
	cout << S.transpose()  << endl;                    // Transpose
	cout << "-------------------------------" << endl;
	cout << "diag(T) = " << endl;
	cout << T.diagonal()   << endl;                    // Only the diagonal of Matrix A
	cout << "-------------------------------" << endl;

	cout << "-------------------------------" << endl;
	cout << "Matrix and vector operations" << endl;
	cout << "-------------------------------" << endl;
	cout << "S*x (built with SparseMatrix) = " << endl;
	SparseVector<double> y(S*x);
	cout << y << endl;
	cout << "-------------------------------" << endl;
	cout << "x^T*T (built with SparseMatrix) = " << endl;
	SparseVector<double> z(x.transpose()*T);
	cout << z << endl;
	cout << "-------------------------------" << endl;
	cout << "S*T (built with SparseMatrix) = " << endl;
	SparseMatrix<double> U(S*T);
	cout << U << endl;
	cout << "-------------------------------" << endl;

	cout << "-------------------------------" << endl;
	cout << "Sparse to dense" << endl;
	cout << "-------------------------------" << endl;
	cout << "B (dense) = S (sparse) = " << endl;
	Matrix<double, 10, 10> B;
	B = MatrixXd(S);
	cout << B << endl;
	cout << "-------------------------------" << endl;
	cout << "S (sparse) = B (dense) = " << endl;
	S = B.sparseView();
	cout << S << endl;
	cout << "-------------------------------" << endl;


	return 0;
}
