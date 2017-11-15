#include"Solveurs.h"
#include <chrono>
#include<fstream>
#include <iostream>

using namespace Eigen;
using namespace std;

int main()
{

  //******************************************************************************
  //                                 initialisation
  //******************************************************************************
  Eigen::MatrixXd A;
  Eigen::MatrixXd M;
  Eigen::MatrixXd B;
  Eigen::MatrixXd Id;
  Eigen::VectorXd b;

  int N;
  int a;


  cout << "choisissez N, la taille de la matrice" << endl;
  cin >> N;


  A.resize(N,N);
  M.resize(N,N);
  B.resize(N,N);
  b.resize(N);

  Id.setIdentity(N,N);


  for(int i(0) ; i<N ; i++)
  {
    b(i)=1;
    for(int j(0) ; j<N ; j++)
    {
      int c = (rand() % 2);
      B(i,j)= c; // pb avec le random il rend toujours le même ordre de 0,1...
    }
  }
  //cout << b << endl;
  cout << B << endl;

  cout << "choisissez a pour constuire A la matrice avec a>0" << endl;
  cin >> a;
  A = a*Id + B.transpose()*B; //construction de A
  cout << A << endl;


  //*****************************************************************************
  //                                    Résolution
  //*****************************************************************************
  double tolerance(0);
  int kmax(0);
  Eigen::VectorXd x;
  Eigen::VectorXd d;
  Eigen::VectorXd x0;
  x.resize(N);
  x0.resize(N);

  for(int i(0) ; i<N ; i++)
  {
    x0(i)=0; // pour avoir une initialisation de x0 comme je ne savais pas quoi mettre
  }

  cout << "iteration max ?" << endl;
  cin >> kmax;

  cout << "Tolérance ?" << endl;
  cin >> tolerance;





  Solveur* test1(0);

  test1 = new Solveur(A,x,b,x0,tolerance,kmax);

  int userChoiceSys; // Choix de l’utilisateur
  cout << "------------------------------------" << endl;
  cout << "Choississez le votre méthode de résolution : " << endl;
  cout << "1) Gradient optimal"<< endl;
  cout << "2) Résidu minimum" << endl;
  cout << "3) Gradient Conjugué" <<endl;
  cout << "4) résidu minimum préconditionné à gauche" << endl;
  //cout << "5) " << endl;
  cin >> userChoiceSys;

  if ((userChoiceSys == 4))
  {
    int userChoiceSys2;
    cout << "------------------------------------" << endl;
    cout << "Choississez votre préconditionneur : " << endl;
    cout << "1) diagonal"<< endl;
    cout << "2) Gauss-Seidel symetrique" << endl;
    cout << "3) Cholesky incomplet" <<endl;
    cin >> userChoiceSys2;
    switch (userChoiceSys2)
    {
      case 1:
      test1->preconddiag(A,M);
      //cout << test.GetM() << endl;
      M=test1->GetM();
      break;

      case 2:
      test1->precondgauss(A,M);
      M=test1->GetM();
      break;

      case 3:
      break;
    }

  }

  switch(userChoiceSys)
  {
    case 1:
    test1->Gradientoptimal(A,x,b,x0,tolerance,kmax);
    cout << test1->Getx()<< endl;
    break;

    case 2:
    test1->Residumini( A,  x,  b, x0 , tolerance, kmax);
    cout << test1->Getx()<< endl;
    break;

    case 3:
    test1->Gradientconjugue( A,  x, b, x0 , tolerance, kmax);
    cout << test1->Getx()<< endl;
    break;

    case 4:
    test1->Residuminiprecond(M,A,x,b,x0,tolerance,kmax);
    cout << test1->Getx() << endl;
    break;
  }


  //******************************************************************************
  //                             tracé de la norme du résidu
  //******************************************************************************









}
