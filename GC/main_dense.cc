#include <iostream>
#include<complex>
#include "Dense"
using namespace Eigen;
using namespace std;

int main()
{
	Matrix<double, 3, 3> A, R;              // Fixed rows and cols. Same as Matrix3d.
	Matrix<double, 3, Dynamic> B;           // Fixed rows, dynamic cols.
	Matrix<double, Dynamic, Dynamic> C;     // Full dynamic. Same as MatrixXd.
	MatrixXcd D;                            // Full dynamic but complex matrix
	Vector3d x, y;                          // 3x1 double matrix.
	VectorXd v;                             // Dynamic column vector of doubles

	// Matrix resize
	// A.resize(4, 4);                      // Runtime error if assertions are on.
	// B.resize(4, 9);                      // Runtime error if assertions are on.
	A.resize(3, 3);                         // Ok; size didn't change.
	B.resize(3, 9);                         // Ok; only dynamic cols changed.
	C.resize(5, 7);                          // Ok; C is dymanic (MatrixXd)
	D.resize(3, 3);


	A << 1, 2, 3,                           // Initialize A. The elements can also be
	4, 5, 6,                                // matrices.
	7, 8, 9;
	B << A, A, A;                           // B is three horizontally stacked A's.
	D.real() = A;
	D.imag() << 10, 11, 12,
	13, 14, 15,
	16, 17, 18;
	x << 3, 8, 9;
	y << 2, 8, 6;

	cout << "-------------------------------" << endl;
	cout << "A = (square and double)" << endl;
	cout << A << endl;                       // Print Matrix A
	cout << "-------------------------------" << endl;
	cout << "B = (double)" << endl;
	cout << B << endl;                       // Print Matrix B
	cout << "-------------------------------" << endl;
	cout << "D = (square and complex)" << endl;
	cout << D << endl;                       // Print Matrix D
	cout << "-------------------------------" << endl;
	cout << "x = (vector)" << endl;
	cout << x << endl;                       // Print Vector x
	cout << "-------------------------------" << endl;
	cout << "y = (vector)" << endl;
	cout << y << endl;                       // Print Vector y

	// Basic usage
	cout << "Rk : The numbering starts at 0" << endl;
	cout << "-------------------------------" << endl;
	int size_x = x.size();                  // size of Vector x
	int i = floor(size_x/2);
	double x_i = x(i);                      // x(i), i = 0...x.size()
	cout << "x(" << i <<") = " << x_i << endl;
	cout << "-------------------------------" << endl;
	int num_rows_B = B.rows();              // number of rows of Matrix B
	int num_cols_B = B.cols();              // number of columns of Matrix B
	i = floor(num_rows_B/2);
	int j = floor(num_cols_B/2);
	double B_i_j = B(i,j);                  // B(i,j), i=0..B.rows()-1, j=0..B.cols()-1
	cout << "B(" << i << ", " << j <<") = " << B_i_j << endl;
	cout << "-------------------------------" << endl;

	// Define a Matrix
	C = MatrixXd::Identity(5,5);            // Identity matrix
	// C = MatrixXd::Zero  (rows,cols);     // Matrix with zeros
	// C = MatrixXd::Ones  (rows,cols);     // Matrix with ones
	// C = MatrixXd::Random(rows,cols);     // MatrixXd::Random returns uniform random
                                          // numbers in (-1, 1).
	cout << "C = (5x5 identity)" << endl;
	cout << C << endl;                      // Print Matrix C
	cout << "-------------------------------" << endl;

	// Set a Matrix
	R.setRandom(R.rows(),R.cols());         // Random matrix
	// A.setZero    (A.rows(),A.cols());    // Matrix with zeros
	// A.setOnes    (A.rows(),A.cols());    // Matrix with ones
	// A.setIdentity(A.rows(),A.cols());    // Identity matrix
	cout << "R = (5x5 random)" << endl;
	cout << R << endl;                      // Print Matrix R
	cout << "-------------------------------" << endl;

	// Linearly Spaced Vector
	v = VectorXd::LinSpaced(11,0.,1.);      // vector = [0. 0.1 ... 1] (size = 11)
	cout << "v = (Linearly Spaced Vector)" << endl;
	cout << v << endl;
	cout << "-------------------------------" << endl;

	// Matrix slicing and blocks.
	cout << "Slicing and blocks" << endl;
	cout << "-------------------------------" << endl;
	int k(3), N(5);
	cout << "v(i) (i = " << k <<" ... " << k+N-1 <<") = " << endl;
	cout << v.segment(k, N) << endl;       // v(k) ... v(k+n-1)
	cout << "-------------------------------" << endl;

	int l(0), m(1), P(3), Q(2);
	cout << "R(" << l << ", " << m << ") ... R(" << l << ", " << m+Q-1 <<")" << endl;
	cout << "   ." << "           ." << endl;
	cout << "   ." << "           ." << endl;
	cout << "   ." << "           ." << endl;
	cout << "R(" << l+P-1 << ", " << m << ") ... R(" << l+P-1 << ", " << m+Q-1 <<")" << endl;
	cout << "=" << endl;
	cout << R.block(l, m, P, Q) << endl;
	cout << "-------------------------------" << endl;

	cout << "Adjoint, transpose, conjugate " << endl;
	cout << "-------------------------------" << endl;
	cout << "A^T = " << endl;
	cout << A.transpose()  << endl;                    // Transpose
	cout << "-------------------------------" << endl;
	cout << "diag(A) = " << endl;
	cout << A.diagonal()   << endl;                    // Only the diagonal of Matrix A
	cout << "-------------------------------" << endl;
	cout << "bar(D)" << endl;
	cout << D.conjugate()  << endl;                    // Conjugate
  cout << "-------------------------------" << endl;
	cout << "D^* = " << endl;
	cout << D.adjoint()    << endl;                    // Adjoint
	cout << "-------------------------------" << endl;

	cout << "Matrix and vector operations" << endl;
	cout << "Operator overloading" << endl;
	cout << "-------------------------------" << endl;
	cout << "A*x = " << endl;
	cout << A*x << endl;
	cout << "-------------------------------" << endl;
	cout << "x^T*A = " << endl;
	cout << x.transpose()*A << endl;
	cout << "-------------------------------" << endl;
	cout << "R*A" << endl;
	cout << R*A << endl;
	cout << "-------------------------------" << endl;

	cout << "Vectorized operations on each element independently" << endl;
	cout << "-------------------------------" << endl;
	cout << "R.*A" << endl;
	cout << R.array()*A.array() << endl;
	auto F = R.array()*A.array() > 0.;
	cout << "-------------------------------" << endl;
	cout << "D = (R.*A > 0) " << endl;
	cout << "  = 1 if R(i)*A(i) > 0" << endl;
	cout <<  F << endl;
	cout << "-------------------------------" << endl;
	cout << "A.^2 != A^2" << endl;
	cout << A.array().square() << endl;
	cout << "different to" << endl;
	cout << A*A << endl;
	cout << "-------------------------------" << endl;
	cout << "cos(A(i))" << endl;
	cout << A.array().cos() << endl;
	cout << "-------------------------------" << endl;

	cout << "Dot products, norms, ... " << endl;
	cout << "norm(x) = " << endl;
	cout << x.norm() << endl;
	cout << "-------------------------------" << endl;
	cout << "x . y = " << endl;
	cout << x.dot(y) << endl;
	cout << "-------------------------------" << endl;
	cout << " x ^ y = " << endl;
	cout << x.cross(y) << endl;
	cout << "-------------------------------" << endl;

	return 0;
}
