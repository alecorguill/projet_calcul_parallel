#include "parallel_conduction.hpp"
#include "util.hpp"

#include <assert.h>
#include <stdio.h>
#include <iostream>

using namespace std;


void test_parse_file(){
  printf("TEST PARSE_FILE...");
  config_t c;
  parse_file("test/config.cfg", c);
  //cout << c.Nx << " " << c.Ny << " " << c.Lx << " " << c.Ly << " " << c.D << " " << c.choix;
  assert(c.Nx == 1000);
  assert(c.Ny == 1000);
  assert(c.Lx == 1.0);
  assert(c.Ly == 1.0);
  assert(c.D == 1);
  assert(c.choix == 1);
  printf("OK\n");
}

int main(){
  test_parse_file();
}
