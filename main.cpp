#include <iostream>
#include <vector>
#include <cstdlib>

#include "./include/system.hpp"

int brute_force(int L) {
  int seg = 2 * L * (L + 1);

  int cnt = 0;

  for(int S = 0; S < (1 << seg); ++S) {
    Configuration Ctmp(seg);
    for(int i = 0; i < seg; ++i) {
      if((S >> i) & 1) Ctmp[i] = true;
    }
    int etmp = energy(Ctmp, L);
    std::cerr << etmp << std::endl;
    if(etmp == 0) ++cnt;
  }

  return cnt;
}


int main(int argc, char* argv[]) {
  std::string program = argv[0];
  if(argc != 2) {
    std::cerr << "Usage: " << program << " L" << std::endl;
  }

  int L = std::atoi(argv[1]);

  int cnt = brute_force(L);

  std::cout << cnt << std::endl;
  return 0;
}
