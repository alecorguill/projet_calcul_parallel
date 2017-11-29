#ifndef PARALLEL_CONDUCTION
#define PARALLEL_CONDUCTION

#include "util.hpp"
#include "Dense"

/* flux */
double f(int i, int j, config_t& c);

/* conditions de bord */
double g(int i, int j, config_t& c);
double h(int i, int j, config_t& c);

/* traduction des indices matricielles en indice linéaires */
void indice(int k, int& i, int& j, config_t& c);

/* renvoie le morceau de domaine en fonction du proc  */
void charge(int me, int &i0, int &i1, config_t& c);

Eigen::VectorXd second_membre(int me, Eigen::VectorXd u, config_t& c);

void Gradientconjugue(Eigen::MatrixXd A, Eigen::VectorXd u, Eigen::VectorXd b,Eigen::VectorXd x0 ,double tolerance, int kmax);

void Remplissage(Eigen::MatrixXd& A, int Nx, int Ny, config_t& c);

#endif
