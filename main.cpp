#include <iostream>
#include <vector>
#include <cstdlib>

#include "./include/system.hpp"
#include "./include/MCMC.hpp"

int brute_force(int L) {
  int seg = 2 * L * (L + 1);

  int cnt = 0;

  for(int S = 0; S < (1 << seg); ++S) {
    if(!(S & 1) || ((S >> (seg / 2) & 1))) continue;
    Configuration Ctmp(seg);
    for(int i = 0; i < seg; ++i) {
      if((S >> i) & 1) Ctmp[i] = true;
    }
    int etmp = energy(Ctmp, L);
    std::cerr << etmp << std::endl;
    if(etmp == 0) ++cnt;
  }

  cnt *= 2;

  return cnt;
}


int main(int argc, char* argv[]) {
  std::string program = argv[0];
  if(argc != 2) {
    std::cerr << "Usage: " << program << " L" << std::endl;
  }

  int L = std::atoi(argv[1]);

  std::vector<double> entropy;
  Wang_Landau(L, entropy);

  double P0_tmp = multicanonical_sampling(L, entropy);
  std::cerr << P0_tmp << std::endl;
  double estimate = P0_tmp / softmax(entropy, 0);
  std::cout << estimate << std::endl;
  return 0;
}
