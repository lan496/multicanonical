#include <iostream>
#include <vector>
#include <cmath>
#include <queue>

#include "../include/system.hpp"

struct Edge {
  int src, dst;
  int weight;
  Edge() { ; }
  Edge(int src, int dst) : src(src), dst(dst), weight(1) { ; }
};

typedef std::vector<Edge> Edges;
typedef std::vector<Edges> Graph;

void convert_configuration_to_graph(const Configuration& C, int L, Graph& G) {
  int V = (L + 1) * (L + 1);
  G.assign(V, Edges());

  for(int i = 0; i < L * (L + 1); ++i) {
    if(!C[i]) continue;
    int src = (i / L) * (L + 1) + i % L;
    int dst = src + 1;
    G[src].push_back(Edge(src, dst));
    G[dst].push_back(Edge(dst, src));
  }
  for(int i = L * (L + 1); i < 2 * L * (L + 1); ++i) {
    if(!C[i]) continue;
    int ii = i - L * (L + 1);
    int src = (ii / (L + 1)) * (L + 1) + ii % (L + 1);
    int dst = src + L + 1;
    G[src].push_back(Edge(src, dst));
    G[dst].push_back(Edge(dst, src));
  }
}

//http://www.prefield.com/algorithm/container/union_find.html
struct UnionFind {
  std::vector<int> data;
  UnionFind(int size) : data(size, -1) { }
  bool unionSet(int x, int y) {
    x = root(x); y = root(y);
    if (x != y) {
      if (data[y] < data[x]) std::swap(x, y);
      data[x] += data[y]; data[y] = x;
    }
    return x != y;
  }
  bool findSet(int x, int y) {
    return root(x) == root(y);
  }
  int root(int x) {
    return data[x] < 0 ? x : data[x] = root(data[x]);
  }
  int size(int x) {
    return -data[root(x)];
  }
};

int energy(const Configuration& C, int L) {
  Graph G;
  convert_configuration_to_graph(C, L, G);

  int V = (L + 1) * (L + 1);
  int e = 0;

  for(int i = 0; i < V; ++i) {
    int degree = G[i].size();
    if(i == 0 || i == V - 1) {
      e += std::abs(degree - 1);
    }else{
      e += std::min(std::abs(degree - 2), degree);
    }
  }

  UnionFind uf(V);
  for(int i = 0; i < V; ++i) {
    for(Edge e: G[i]) {
      uf.unionSet(e.src, e.dst);
    }
  }
  for(int i = 0; i < V; ++i) {
    if(G[i].size() > 0 && !uf.findSet(0, i)) ++e;
  }

  return e;
}

int proposed_energy(const Configuration& C, int L, int pos) {
  Configuration Ctmp(C);
  Ctmp[pos] = !Ctmp[pos];
  int e = energy(Ctmp, L);
  return e;
}
