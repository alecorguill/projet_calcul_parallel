#ifndef PARALLEL_CONDUCTION
#define PARALLEL_CONDUCTION

#include "util.hpp"

/* flux */
double f(int i, int j, config_t& c);

/* conditions de bord */
double g(int i, int j, config_t& c);
double h(int i, int j, config_t& c);

/* traduction des indices matricielles en indice linéaires */
void indice(int k, int& i, int& j, config_t& c);

/* renvoie le morceau de domaine en fonction du proc  */
void charge(int me, int &i0, int &i1, config_t& c);

#endif
