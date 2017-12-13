#include <iostream>
#include <vector>
#include <cmath>
#include <random>

#include "../include/system.hpp"

bool is_flat(const std::vector<int>& hist) {
  int N = hist.size();

  bool ret = true;
  int cnt = 0;
  for(int i = 0; i < N - 1; ++i) {
    if(!hist[i] || !hist[i + 1]) continue;
    if((double)std::abs(hist[i] - hist[i + 1]) / hist[i] > 0.1) ret = false;
    ++cnt;
  }

  if(cnt < N / 3) ret = false;
  return ret;
}

void Wang_Landau(int L, int delta, std::vector<double> entropy) {
  int V = (L + 1) * (L + 1);

  int Emax = 3 * V;
  int Emin = 0;

  int N = (Emax - Emin + delta - 1) / delta;
  std::vector<int> hist(N, 0);
  entropy.assign(N, 0);
  double f = 1.0;

  double eps = 1e-4;

  std::random_device rnd;
  std::mt19937 mt(rnd());
  std::uniform_real_distribution<double> bern(0.0, 1.0);

  Configuration current(V);
  for(int i = 0; i < V; ++i) {
    current[i] = mt() % 2;
  }
  int current_e = energy(current, L);

  while(f > eps) {
    int propose_idx = mt() % V;
    int proposed_e = proposed_energy(current, L, propose_idx);

    if(bern(mt) < std::exp(entropy[current_e] - entropy[proposed_e])) {
      current_e = proposed_e;
      current[propose_idx] = !current[propose_idx];
    }

    hist[current_e / delta] += 1;
    entropy[current_e / delta] += f;

    for(int i = 0; i < V; ++i) std::cerr << " " << entropy[i];
    std::cerr << std::endl;

    if(is_flat(hist)) {
      for(int i = 0; i < N; ++i) {
        hist[i] = 0;
      }
      f *= 0.5;
    }
  }
}
