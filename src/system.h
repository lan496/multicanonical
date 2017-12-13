#ifndef SYSTEM_H
#define SYSTEM_H

#include <iostream>
#include <vector>

class Edge {
  int src, dst;
  Graph() : {}
  Graph(src, dst) : src(src), dst(dst) {}
};

typedef std::vector<std::vector<Edge>> Graph;
typedef std::vector<bool> Configuration;

void convert_configuration_to_graph(const Configuration&, Graph&);
double energy(const Configuration&);

#endif
