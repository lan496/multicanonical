#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>

#include "./include/system.hpp"
#include "./include/MCMC.hpp"

void brute_force(int L, std::vector<int>& hist) {
  int seg = 2 * L * (L + 1);

  int V = (L + 1) * (L + 1);
  hist.assign(3 * V, 0);

  for(int S = 0; S < (1 << seg); ++S) {
    Configuration Ctmp(seg);
    for(int i = 0; i < seg; ++i) {
      if((S >> i) & 1) Ctmp[i] = true;
    }
    int etmp = energy(Ctmp, L);
    ++hist[etmp];
    std::cerr << etmp << std::endl;
  }
}


int main(int argc, char* argv[]) {
  std::string program = argv[0];
  if(argc != 2) {
    std::cerr << "Usage: " << program << " L" << std::endl;
  }

  int L = std::atoi(argv[1]);

  std::vector<double> entropy;
  Wang_Landau(L, entropy);

  int t = 1000;
  long double mean = 0, mean2 = 0;

  for(int i = 0; i < t; ++i) {
    long double tmp = multicanonical_sampling(L, entropy);
    std::cerr << "  " << tmp << std::endl;

    if(i == 0) {
      mean = tmp;
      mean2 = tmp * tmp;
    }else{
      mean = (1.0 - 1.0 / i) * mean + tmp / i;
      mean2 = (1.0 - 1.0 / i) * mean2 + tmp * (tmp / i);
      double err = std::sqrt((mean2 - mean * mean) / i);
      std::cout << mean << " " << err << " " << err / mean << std::endl;
    }

  }
  return 0;
}
