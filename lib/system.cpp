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

int dijkstra(const Graph& G, int s, int t) {
  int V = G.size();
  int INF = V;
  std::vector<int> d(V, INF);
  d[s] = 0;

  using P = std::pair<int,int>;
  std::priority_queue<P, std::vector<P>, std::greater<P>> que;
  que.push(P(0, s));
  while(!que.empty()) {
    P q = que.top();
    que.pop();
    if(d[q.second] < q.first) continue;
    for(Edge e: G[q.second]) {
      if(d[e.dst] > d[e.src] + e.weight) {
        d[e.dst] = d[e.src] + e.weight;
        que.push(P(d[e.dst], e.dst));
      }
    }
  }

  return d[t];
}

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

  int E = 0;
  for(int i = 0; i < V; ++i) {
    E += G[i].size();
  }
  E /= 2;

  int dist = dijkstra(G, 0, V - 1);
  e += std::abs(dist - E);

  return e;
}

int proposed_energy(const Configuration& C, int L, int pos) {
  Configuration Ctmp(C);
  Ctmp[pos] = !Ctmp[pos];
  int e = energy(Ctmp, L);
  return e;
}
