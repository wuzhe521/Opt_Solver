#include "../core/matrix.h"
#include "../core/vector.h"

using namespace gons::utilites::LOG_MSG;

int main() {
  gons::Matrix<double, 2, 2> matrix = {{1, 2}, {3, 4}};
  LOG(matrix);
  gons::Vector<double, 2> vector = {1, 2};
  LOG(vector);
  
  vector = matrix * vector.Transpose();
  LOG(vector);

  LOG(vector.Norm2());
  return 0;
}