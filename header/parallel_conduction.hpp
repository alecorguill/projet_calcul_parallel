#ifndef PARALLEL_CONDUCTION
#define PARALLEL_CONDUCTION

#include "util.hpp"
#include "Eigen"

/* flux */
double f(int i, int j, config_t& c);

/* conditions de bord */
double g(int i, int j, config_t& c);
double h(int i, int j, config_t& c);

/* traduction des indices matricielles en indice linéaires */
void indice(int k, int& i, int& j, config_t& c);

/* renvoie le morceau de domaine en fonction du proc  */
void charge(int me, int &i0, int &i1, config_t& c);

// Fonction permettant la modification du second membre pour envoyer les informations de recouvrement
// Elle prend en paramètre le processeur et le vecteur u et renvoie un vecteur local du second membre
Eigen::VectorXd second_membre(int me, Eigen::VectorXd u, config_t& c);

// Fonction permettant de vérifier la Convergence pour sortir des itérations de Schwarz
// Elle prend en paramètre 2 vecteurs et une tolérance et renvoie un nombre
// Si le processeur a une bonne convergence on renvoie 1 sinon c'est 0
int Convergence(Eigen::VectorXd utemp, Eigen::VectorXd u, double e);

// Procédure pour faire le gradient conjugué
// Prend en paramètre une matrice A, un vecteur u et un second membre b : on résoud Au = b
// Le vecteur x0 est nécessaire pour le gradient conjugué, dans notre cas il sera prit égale à 0
// Tolerance et kmax sont la pour gérer la sortie du gradient conjugué, on tend vers une tolérance et un garde fou kmax qui est un nombre d'itérations maximas
// La résolution se fait en locale, nécessitant ainsi la présence du me
void Gradientconjugue(Eigen::MatrixXd A, Eigen::VectorXd& u, Eigen::VectorXd b,Eigen::VectorXd x0 ,double tolerance, int kmax, config_t& c, int me);

// Fonction Bi-Gradient Conjugué permettant de faire comme le gradient conjugué pour des matrices non symétriques
void BIGradientconjugue(Eigen::MatrixXd A, Eigen::VectorXd& u, Eigen::VectorXd b,Eigen::VectorXd x0 ,double tolerance, int kmax, config_t& c, int me);

// Procédure permettant de remplir la matrice A connaissant son Nx et son Ny 
void Remplissage(Eigen::MatrixXd& A, int Nx, int Ny, config_t& c);

#endif
