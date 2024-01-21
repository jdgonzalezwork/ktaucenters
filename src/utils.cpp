#include <Rcpp.h>
#include <algorithm>
#include <queue>
using namespace Rcpp;

double median_cpp(NumericVector x) {

  //' Computes the median of a numeric vector
  //'
  //' @param x numeric vector.
  //'
  //' @return median of a numeric vector.
  std::size_t size = x.size();
  std::sort(x.begin(), x.end());
  if (size % 2 == 0)
    return (x[size / 2 - 1] + x[size / 2]) / 2.0;
  return x[size / 2];
}

class IndexComparator {
public:
  IndexComparator(NumericVector data_, bool decreasing_)
      : data(data_.begin()), decreasing(decreasing_) {}

  inline bool operator()(int i, int j) const {
    if (decreasing) {
      // Ties: If there are ties, then the first value is returned and the other
      // one is ignored.
      return data[i] > data[j] || (data[i] == data[j] && j > i);
    } else {
      return data[i] < data[j] || (data[i] == data[j] && j > i);
    }
  }

private:
  NumericVector::iterator data;
  bool decreasing;
};

class IndexQueue {
public:
  typedef std::priority_queue<int, std::vector<int>, IndexComparator> Queue;

  IndexQueue(NumericVector data_, bool decreasing_)
      : comparator(data_, decreasing_), q(comparator), data(data_) {}

  inline operator IntegerVector() {
    int n = q.size();
    IntegerVector res(n);
    for (int i = n - 1; i >= 0; --i) {
      // Add +1 for 1-based R indexing
      res[i] = q.top();
      q.pop();
    }
    return res;
  }
  inline void input(int i) {
    // if( data[ q.top() ] < data[i] ) {
    if (comparator(i, q.top())) {
      q.pop();
      q.push(i);
    }
  }
  inline void pop() { q.pop(); }
  inline void push(int i) { q.push(i); }

private:
  IndexComparator comparator;
  Queue q;
  NumericVector data;
};

IntegerVector top_index(NumericVector v, int n, bool decreasing) {

  //' Find the first nth indices from a vector
  //'
  //' @param v a numeric vector from which indices are found.
  //' @param n number of indices to find.
  //' @param decreasing boolean value indicating if the vector must be in
  //' decreasing order
  //'
  //' @return
  //' A vector with first nth indices.
  int size = v.size();

  // not interesting case. Less data than n
  if (size < n) {
    return seq(0, n - 1);
  }

  IndexQueue q(v, decreasing);
  for (int i = 0; i < n; ++i)
    q.push(i);
  for (int i = n; i < size; ++i)
    q.input(i);
  return q;
}

double max_tolerance(NumericMatrix x, NumericMatrix y) {

  //' Estimates maximum tolerance between two given matrix row wise
  //'
  //' @param x numeric matrix.
  //' @param y numeric matrix.
  //'
  //' @return max euclidean difference two given matrix row wise.

  if ((x.nrow() != y.nrow()) || (x.ncol() != y.ncol())) {
    stop("Both matrices must have same dimension");
  }
  const std::size_t n = x.nrow();
  NumericVector tolerances(no_init(n));

  for (std::size_t i = 0; i < n; ++i) {
    tolerances[i] = sqrt(sum(pow(x.row(i) - y.row(i), 2)));
  }

  return max(tolerances);
}
