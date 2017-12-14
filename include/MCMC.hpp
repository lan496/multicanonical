#ifndef MCMC_H
#define MCMC_H

#include <vector>

void Wang_Landau(int, std::vector<double>&);
double multicanonical_sampling(int, const std::vector<double>&);
double softmax(const std::vector<double>&, int);

#endif
