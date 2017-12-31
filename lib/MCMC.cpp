#include <iostream>
#include <vector>
#include <cmath>
#include <random>

#include "../include/system.hpp"

double softmax(const std::vector<double>& entropy, int pos) {
  int nonzero = 0;
  double ave = 0;
  for(double s: entropy) {
    if(s > 0) {
      ++nonzero;
      ave += s;
    }
  }
  ave /= nonzero;

  double Z = 0;
  for(double s: entropy) {
    if(s == 0) continue;
    Z += std::exp(s - ave);
  }

  double ret;
  if(entropy[pos] == 0) {
    ret = 0;
  }else{
    ret = std::exp(entropy[pos] - ave) / Z;
  }

  return ret;
}

bool is_flat(const std::vector<int>& hist) {
  int N = hist.size();

  bool ret = true;
  int cnt = 0;
  for(int i = 0; i < N - 1; ++i) {
    if(!hist[i] || !hist[i + 1]) continue;
    if((double)std::abs(hist[i] - hist[i + 1]) / hist[i] > 1e-2) ret = false;
    ++cnt;
  }

  if(cnt < N / 3) ret = false;
  return ret;
}

void Wang_Landau(int L, std::vector<double>& entropy) {
  int V = (L + 1) * (L + 1);

  int Emax = 2 * V;
  int Emin = 0;

  int N = Emax - Emin;
  std::vector<int> hist(N, 0);
  entropy.assign(N, 0);
  double f = 1.0;

  double eps = 1e-4;

  std::random_device rnd;
  std::mt19937 mt(rnd());
  std::uniform_real_distribution<double> bern(0.0, 1.0);

  int seg = 2 * L * (L + 1);
  Configuration current(seg);
  for(int i = 0; i < seg; ++i) {
    current[i] = mt() % 2;
  }
  int current_e = energy(current, L);

  while(f > eps) {
    int propose_idx = mt() % seg;
    int proposed_e = proposed_energy(current, L, propose_idx);

    if(bern(mt) < std::exp(entropy[current_e] - entropy[proposed_e])) {
      current_e = proposed_e;
      current[propose_idx] = !current[propose_idx];
    }

    hist[current_e] += 1;
    entropy[current_e] += f;

    if(is_flat(hist)) {
      for(double e : entropy) std::cerr << " " << e;
      std::cerr << std::endl;
      std::cerr << f << std::endl;

      for(int i = 0; i < N; ++i) {
        hist[i] = 0;
      }
      f *= 0.5;
    }
  }
}

long double estimate_number(const std::vector<double>& entropy, const std::vector<long int>& hist, int seg) {
  int nonzero = 0;
  long double ave = 0;
  for(double s: entropy) {
    if(s > 0) {
      ++nonzero;
      ave += s;
    }
  }
  ave /= nonzero;

  long double total = 0.0;
  for(int i = 0; i < (int)entropy.size(); ++i) {
    if(entropy[i] == 0) continue;
    total += hist[i] * std::exp(entropy[i] - ave);
  }

  long double ret = hist[0] * std::exp(entropy[0] - ave) / total;
  for(int i = 0; i < seg; ++i) {
    ret *= 2;
  }

  return ret;
}

long double multicanonical_sampling(int L, const std::vector<double>& entropy) {
  std::random_device rnd;
  std::mt19937 mt(rnd());
  std::uniform_real_distribution<double> bern(0.0, 1.0);

  int seg = 2 * L * (L + 1);
  Configuration current(seg);
  for(int i = 0; i < seg; ++i) {
    current[i] = mt() % 2;
  }
  int current_e = energy(current, L);

  int n_sample = 100000;
  int mcs = 4;

  int V = (L + 1) * (L + 1);
  std::vector<long int> hist(2 * V, 0);

  for(int i = 0; i < n_sample / 100; ++i) {
    int propose_idx = mt() % seg;
    int proposed_e = proposed_energy(current, L, propose_idx);

    if(bern(mt) < std::exp(entropy[current_e] - entropy[proposed_e])) {
      current_e = proposed_e;
      current[propose_idx] = !current[propose_idx];
    }
  }

  for(int i = 0; i < mcs * n_sample; ++i) {
    int propose_idx = mt() % seg;
    int proposed_e = proposed_energy(current, L, propose_idx);

    if(bern(mt) < std::exp(entropy[current_e] - entropy[proposed_e])) {
      current_e = proposed_e;
      current[propose_idx] = !current[propose_idx];
    }

    if(i % mcs == 0) ++hist[current_e];
  }

  long double estimate = estimate_number(entropy, hist, seg);
  return estimate;
}

