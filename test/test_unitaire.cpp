#include <parallel_conduction.hpp>
#include <assert.h>
#include <stdio.h>

void test_indice(){
  printf("TEST INDICE...");
  int nx=5,ny=5;
  int i,j,k;
  k = 11;  
  indice(k,nx,ny,i,j);
  assert(i == 2); assert(j == 1);
  k = 24;
  indice(k,nx,ny,i,j);
  assert(i == 4); assert(j == 4);
  printf("OK\n");
}

void test_charge(){
  printf("TEST CHARGE...");
  int ny=17;
  int me,np, i0, i1;
  np = 3;
  charge(0,ny,np,i0,i1);
  assert(i0 == 0); assert(i1 == 5);
  charge(1,ny,np,i0,i1);
  assert(i0 == 6); assert(i1 == 11);
  charge(2,ny,np,i0,i1);
  assert(i0 == 12); assert(i1 == 16);
  printf("OK\n");
}
int main(){
  test_indice();
  test_charge();
}