#include <iostream>
#include <fstream>
#include <limits>
#include <queue>
#include <unordered_map>
#include <vector>

template <class Vertex>
struct Edge {
  Edge(const Vertex& from, const Vertex& to)
      : vertices_(from, to) {}

  const Vertex& GetStart() const { return vertices_.first; }
  const Vertex& GetTarget() const { return vertices_.second; }

 private:
  std::pair<Vertex, Vertex> vertices_;
};

template <typename Vertex = int>
class Graph {
 public:
  using VertexT = Vertex;
  using EdgeT = Edge<VertexT>;

  Graph(const std::vector<EdgeT>& edges) {
    for (const auto& edge : edges) {
      adjacent_[edge.GetStart()].emplace_back(edge);
    }
  }

  const std::vector<EdgeT>& GetOutgoingEdges(const Vertex& vertex) {
    return adjacent_[vertex];
  }

 private:
  std::unordered_map<VertexT, std::vector<EdgeT>> adjacent_;
};

template <class Graph, class Visitor>
void BFS(Graph& graph, const typename Graph::VertexT& start,
              Visitor& visitor) {
  using Vertex = typename Graph::VertexT;
  std::unordered_map<Vertex, int> dist_eval;
  std::queue<std::pair<int, Vertex>> bfs_queue;
  dist_eval[start] = 0;
  bfs_queue.push({dist_eval[start], start});
  while (!bfs_queue.empty()) {
    auto from = bfs_queue.front().second;
    bfs_queue.pop();
    for (const auto& edge : graph.GetOutgoingEdges(from)) {
      const auto& to = edge.GetTarget();
      if (!dist_eval.contains(to) || dist_eval[from] + 1 < dist_eval[to]) {
        dist_eval[to] = dist_eval[from] + 1;
        bfs_queue.push({dist_eval[to], to});
      }
    }
  }
  visitor.Fill(dist_eval);
}

template <class Graph>
class AbstractBFSVisitor {
 public:
  virtual void DiscoverEdge(const Graph::EdgeT& vertex) = 0;
  virtual ~AbstractBFSVisitor() = default;
  std::vector<typename Graph::VertexT> answer;
};

template <class Graph>
class ShortestPathsVisitor : public AbstractBFSVisitor<Graph> {
 public:
  using VertexT = Graph::VertexT;
  using ParrentClass = AbstractBFSVisitor<Graph>;
  using ParrentClass::answer;

  void DiscoverEdge(const Graph::EdgeT& edge) override {
    ancestors_[edge.GetTarget()] = edge.GetStart();
    answer[edge.GetTarget()] = std::min(
        answer[edge.GetTarget()], answer[edge.GetStart()] + 1);
  }
  void CreateAns(int vertexes, int start) {
    answer.resize(vertexes);
    for (int i = 0; i < vertexes; ++i) {
      answer[i] = std::numeric_limits<int>::max();
    }
    answer[start] = 0;
  }

  void Fill(std::unordered_map<VertexT, int>& ans) {
    for (auto& elem : ans) {
      answer[elem.first] = elem.second;
    }
  }

  int GetAns(int ind) { return answer[ind]; }

 private:
  std::unordered_map<VertexT, VertexT> ancestors_;
};

ShortestPathsVisitor<Graph<>> CallBFS(std::vector<Edge<int>>& all_edges,
                                           int vertexes, int start) {
  Graph graph(all_edges);
  ShortestPathsVisitor<Graph<>> visitor;
  visitor.CreateAns(vertexes, start);
  BFS(graph, start, visitor);
  return visitor;
}

void Print(int vertexes, ShortestPathsVisitor<Graph<>>& visitor) {
  for (int elem = 0; elem < vertexes; ++elem) {
    std::cout << visitor.GetAns(elem) << " ";
  }
  std::cout << "\n";
}

void GetData() {
  try {
    std::ifstream file("graph.txt");
    if (!file) {
      throw std::runtime_error("Ошибка открытия файла");
    }
    std::vector<Edge<int>> all_edges;
    int vertexes;
    int edges;
    file >> vertexes >> edges;
    if (file.fail()) {
      throw std::runtime_error("Ошибка чтения количества вершин и рёбер");
    }
    for (int j = 0; j < edges; ++j) {
        int fst;
        int scd;
        int value;
        file >> fst >> scd;
        Edge edge_fst(fst, scd);
        Edge edge_scd(scd, fst);
        if (file.fail()) {
          throw std::runtime_error("Ошибка чтения рёбер графа");
        }
        all_edges.emplace_back(edge_fst);
        all_edges.emplace_back(edge_scd);
    }
    int start;
    file >> start;
    if (file.fail()) {
      throw std::runtime_error("Ошибка чтения начальной вершины");
    }
    ShortestPathsVisitor<Graph<>> visitor =
        CallBFS(all_edges, vertexes, start);
    Print(vertexes, visitor);
  }
  catch (const std::runtime_error& e) {
    std::cerr << "Ошибка: " << e.what() << "\n";
  }
}

int main() { GetData(); }