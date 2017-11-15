#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>

#include "Dense"
#include "Sparse"

class Solveur
{
  private: // les attributs de la classe
    //le système à résoudre est A*x=b
    Eigen::MatrixXd _A; //matrice du système à résoudre
    Eigen::MatrixXd _M; //matrice du préconditinneur
    Eigen::VectorXd _x;  //solution du système à résoudre
    Eigen::VectorXd _b;  // second membre du système à résoudre
    Eigen::VectorXd _x0; // l'opérateur peut choisir l'initialisation

    double _tolerance;   // l'opérateur peut choisir la tolérance qu'il souhaite
    int _kmax;           // l'opérateur peut choisir le nombre d'opération maximale


  public: //
  //constructeur
  Solveur(Eigen::MatrixXd A, Eigen::VectorXd x, Eigen::VectorXd b,Eigen::VectorXd x0 ,double tolerance, int kmax);


  //méthode du gradient à pas optimal
  void Gradientoptimal(Eigen::MatrixXd A, Eigen::VectorXd x, Eigen::VectorXd b,Eigen::VectorXd x0 ,double tolerance, int kmax);
    //méthode du résidu minimum
  void Residumini(Eigen::MatrixXd A, Eigen::VectorXd x, Eigen::VectorXd b,Eigen::VectorXd x0 ,double tolerance, int kmax);
  //méthode du gradient conjugué
  void Gradientconjugue(Eigen::MatrixXd A, Eigen::VectorXd x, Eigen::VectorXd b,Eigen::VectorXd x0 ,double tolerance, int kmax);

  //méthode du résidu minimum préconditionné à gauche
  void Residuminiprecond(Eigen::MatrixXd M,Eigen::MatrixXd A, Eigen::VectorXd x, Eigen::VectorXd b,Eigen::VectorXd x0 ,double tolerance, int kmax);

  //pour avoir _x
  const Eigen::VectorXd & Getx() const;
//pour avoir _x
  const Eigen::MatrixXd & GetM() const;
  //préconditionneur diag
  void preconddiag(Eigen::MatrixXd A,Eigen::MatrixXd M);
  //préconditionneur Gauss-Seidel
  void precondgauss(Eigen::MatrixXd A,Eigen::MatrixXd M);
  //préconditionneur Cholesky incomplet
  /*void precondchol(Eigen::MatrixXd A,Eigen::MatrixXd M);
*/

  //Pour avoir la norme du résidu
  //double Getrk();
};
