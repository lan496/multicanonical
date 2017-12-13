#ifndef SYSTEM_H
#define SYSTEM_H

#include <iostream>
#include <vector>

typedef std::vector<bool> Configuration;

int energy(const Configuration&, int L);
int proposed_energy(const Configuration&, int L, int pos);

#endif
