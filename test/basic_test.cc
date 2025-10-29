#include "../core/matrix.h"
#include "../core/utilites.h"
#include "../core/vector.h"


#include <iostream>
using namespace gons::utilites::LOG_MSG;
int main() {

  LOG("GONS Matrix Test Program. Version: " << gons::utilites::GetVersion());
  gons::Matrix<int, 3, 3> mat;
  gons::Matrix<int, 4, 4> mat_b; // Copy constructor
  mat(0, 0) = 1;
  mat(1, 1) = 2;
  mat(2, 2) = 3;
  mat.Print();
  mat(0, 0) = 4;
  mat(1, 1) = 5;
  mat(2, 2) = 6;
  mat.Print("Matrix contents:\n");
  SHW(mat(0, 0));
  // cpoy constructor test
  gons::Matrix<int, 3, 3> mat_copy(mat);
  mat_copy.Print("Copied Matrix contents:\n");
  // INITIALIZER LIST TEST
  gons::Matrix<int, 2, 2> mat_init{{1, 2}, {3, 4}};
  mat_init.Print("Initializer List Matrix contents:\n");

  mat_b.Print("Matrix B contents after assignment:\n");

  std::cout << "Element at (1,2): " << mat(1, 2) << std::endl;

  gons::Matrix<float, 3, 3> A{{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};
  gons::Matrix<float, 3, 3> B{{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};
  gons::Matrix<float, 3, 3> C = A + B;
  C.Print("A + B = \n");

  gons::Matrix<float, 3, 3> D = A - B;
  D.Print("A - B = \n");
  // Matrix multiplication test
  gons::Matrix<float, 2, 3> M1{{1, 2, 3}, {4, 5, 6}};
  gons::Matrix<float, 3, 2> M2{{7, 8}, {9, 10}, {11, 12}};
  gons::Matrix<float, 2, 2> M3 = M1 * M2;
  M3.Print("M1 * M2 = \n");
  // error matrix multiplication test
  gons::Matrix<float, 4, 4> M4 = {
      {1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16}};
  // auto M5 = M1 * M4; // This should trigger a CHECK error
  // transpose test

  // multiply scalar

  gons::Matrix<float, 2, 3> M{{1, 2, 3}, {4, 5, 6}};
  M.Print("M : \n");

  gons::Matrix<float, 3, 2> MT = M.Transpose();
  MT.Print("M' : \n");

  auto Identy = M4.eye();
  Identy.Print("M4 eye() : \n");

  auto M5 = gons::Matrix<float, 4, 4>::Identify();
  M5.Print("M5_I Identify() :");
  LOG(" Normal ");
  LOG_WARNING(" Warnning ")
  LOG_ERROR(" ERROR ");

  auto M6 = M5.Augment();
  M6.Print("M5_I Augment() : \n");

  // inverse test
  gons::Matrix<float, 2, 2> M8{{1, 2}, {3, 4}};
  gons::Matrix<float, 2, 2> M9 = M8.Inverse();
  M9.Print("M8 Inverse() : \n");

   // vector multiplication test
  gons::Matrix<float, 2, 2> M10{{1, 2}, {3, 4}};
  gons::Vector<float, 2> v{5, 6};
  v.Print("v : \n");
  gons::Vector<float, 2> v2 = M10 * v;
  v2.Print("M10 * v : \n");
  return 0;
}