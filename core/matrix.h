#ifndef GONS_CORE_MATRIX_H_
#define GONS_CORE_MATRIX_H_

#include "config.h"
#include "utilites.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace gons {
using namespace utilites::LOG_MSG;

template<typename T, GONS_UINT R>
class Vector;

template <typename T, GONS_UINT R, GONS_UINT C> class Matrix {
protected:
  T data_[R][C] = {0};

public:
  // Default constructor
  Matrix() {}
  // Copy constructor
  Matrix(const Matrix &other) {
    for (GONS_UINT i = 0; i < R; ++i) {
      for (GONS_UINT j = 0; j < C; ++j) {
        data_[i][j] = other.data_[i][j];
      }
    }
  }
  // Assignment operator
  Matrix &operator=(const Matrix &other) {
    if (this != &other) {
      for (GONS_UINT i = 0; i < R; ++i) {
        for (GONS_UINT j = 0; j < C; ++j) {
          data_[i][j] = other.data_[i][j];
        }
      }
    }
    return *this;
  }
  // initializer list constructor
  Matrix(const std::initializer_list<std::initializer_list<T>> &init_list) {
    GONS_UINT i = 0;
    for (const auto &row : init_list) {
      GONS_UINT j = 0;
      for (const auto &elem : row) {
        data_[i][j] = elem;
        ++j;
        CHECK(j > C, "Initializer list has more columns than matrix");
      }
      ++i;
      CHECK(i > R, "Initializer list has more rows than matrix");
    }
  }
  Matrix(const std::initializer_list<T> &init_list) {
    GONS_UINT i = 0;
    GONS_UINT j = 0;
    for (const auto &elem : init_list) {
      data_[i][j] = elem;
      ++j;
      if (j == C) {
        ++i;
        j = 0;
      }
      CHECK(i > R, "Initializer list has more rows than matrix");
    }
  }

public:
  //  Get pointer to the underlying data
  T *get() { return &data_[0][0]; }
  // Getters for number of rows and columns
  GONS_UINT rows() const { return R; }
  GONS_UINT cols() const { return C; }
  std::pair<GONS_UINT, GONS_UINT> shape() const { return {R, C}; }

public:
  // OverLoading the () operator for easy access to matrix elements
  T &operator()(GONS_UINT row, GONS_UINT col) { return data_[row][col]; }
  // Const version of the () operator
  const T &operator()(GONS_UINT row, GONS_UINT col) const {
    return data_[row][col];
  }
  // Overloading the + operator for matrix addition
  Matrix operator+(const Matrix &other) const {
    Matrix result;
    for (GONS_UINT i = 0; i < R; ++i) {
      for (GONS_UINT j = 0; j < C; ++j) {
        result(i, j) = (*this)(i, j) + other(i, j);
      }
    }
    return result;
  }
  // Overloading the - operator for matrix subtraction
  Matrix operator-(const Matrix &other) const {
    Matrix result;
    for (GONS_UINT i = 0; i < R; ++i) {
      for (GONS_UINT j = 0; j < C; ++j) {
        result(i, j) = (*this)(i, j) - other(i, j);
      }
    }
    return result;
  }

public:
  template <GONS_UINT R2, GONS_UINT C2>
  Matrix<T, R, C2> operator*(const Matrix<T, R2, C2> &B) const {
    CHECK(this->cols() != B.rows(),
          "Matrix multiplication dimension mismatch: " +
              std::to_string(this->cols()) + " != " + std::to_string(B.rows()));
    Matrix<T, R, C2> result;
    for (GONS_UINT i = 0; i < this->rows(); ++i) {
      for (GONS_UINT j = 0; j < B.cols(); ++j) {
        result(i, j) = 0;
        for (GONS_UINT k = 0; k < C; ++k) {
          result(i, j) += (*this)(i, k) * B(k, j);
        }
      }
    }
    return result;
  }
  // Scalar multiplication with scalar on the right
  Matrix operator*(const T &scalar) const {
    Matrix result;
    for (GONS_UINT i = 0; i < R; ++i) {
      for (GONS_UINT j = 0; j < C; ++j) {
        result(i, j) = (*this)(i, j) * scalar;
      }
    }
    return result;
  }
  // Scalar multiplication with scalar on the left
  friend Matrix<T, R, C> operator*(const T &scalar,
                                   const Matrix<T, R, C> &matrix) {
    return matrix * scalar;
  }
  // Transpose of the matrix
  Matrix<T, C, R> Transpose() const {
    Matrix<T, C, R> result;
    for (GONS_UINT i = 0; i < R; ++i) {
      for (GONS_UINT j = 0; j < C; ++j) {
        result(j, i) = (*this)(i, j);
      }
    }
    return result;
  }

  Matrix<T, R, R> eye() const {
    CHECK(this->rows() != this->cols(), "this is not a square matrix!")
    Matrix<T, R, R> result;
    const GONS_UINT r = this->rows();
    const GONS_UINT c = this->cols();

    for (GONS_UINT i = 0; i < r; ++i) {
      for (GONS_UINT j = 0; j < c; ++j) {
        if (i == j) {
          result(i, j) = T(1);
        }
      }
    }
    return result;
  }
  // Identity Matrix
  static Matrix<T, R, C> Identity() {
    CHECK(R != C, "Not Square!");
    Matrix<T, R, C> result;
    for (GONS_UINT i = 0; i < R; ++i) {
      for (GONS_UINT j = 0; j < C; ++j) {
        if (i == j) {
          result(i, j) = T(1);
        }
      }
    }
    return result;
  }
  // augment matrix
  Matrix<T, R, 2 * C> Augment() const {
    Matrix<T, R, 2 * C> result;
    for (GONS_UINT i = 0; i < R; ++i) {
      for (GONS_UINT j = 0; j < C; ++j) {
        result(i, j) = (*this)(i, j);
      }
      result(i, C + i) = T(1);
    }
    return result;
  }
  // Inverse of the matrix
  // Guass-Jordan elimination
  // Always find the max in the column
  Matrix<T, R, C> Inverse() const {
    CHECK(R == 0 || C == 0, "Matrix is empty!");
    CHECK(R != C, "Not Square!");
    auto findMaxInColunmIndex = [](GONS_UINT col, Matrix<T, R, 2 * C> &mat) {
      T max = std::abs(mat(col, col));
      GONS_UINT max_row = col;
      for (GONS_UINT i = col; i < R; ++i) {
        if (std::abs(mat(i, col)) > max) {
          max = std::abs(mat(i, col));
          max_row = i;
        }
      }
      return max_row;
    };
    auto swapRow = [](GONS_UINT row1, GONS_UINT row2,
                      Matrix<T, R, 2 * C> &mat) {
      for (GONS_UINT i = 0; i < mat.cols(); ++i) {
        T temp = mat(row1, i);
        mat(row1, i) = mat(row2, i);
        mat(row2, i) = temp;
      }
    };
    auto scaleRow = [](GONS_UINT row, T scale, Matrix<T, R, 2 * C> &mat) {
      for (GONS_UINT i = 0; i < mat.cols(); ++i) {
        mat(row, i) *= scale;
      }
    };
    auto addRow = [](GONS_UINT row1, GONS_UINT row2, T scale,
                     Matrix<T, R, 2 * C> &mat) {
      for (GONS_UINT i = 0; i < mat.cols(); ++i) {
        mat(row1, i) += scale * mat(row2, i);
      }
    };

    Matrix<T, R, C> result;
    Matrix<T, R, 2 *C> Augmented = this->Augment();
    for (GONS_UINT i = 0; i < R; ++i) {
      GONS_UINT max_row = findMaxInColunmIndex(i, Augmented);
      if (max_row != i) {
        swapRow(i, max_row, Augmented);
      }
      CHECK(std::abs(Augmented(i, i)) < GONS_FLT_EPSILON, "Singular matrix!");
      T scale = T(1) / Augmented(i, i);
      scaleRow(i, scale, Augmented);
      for (GONS_UINT j = 0; j < R; ++j) {
        if (j == i)
          continue;
        T eleminator = -Augmented(j, i);
        addRow(j, i, eleminator, Augmented);
      }
    }
    for (GONS_UINT i = 0; i < R; ++i) {
      for (GONS_UINT j = R; j < 2 * R; ++j) {
        result(i, j - R) = Augmented(i, j);
      }
    }
    return result;
  }

  // Frobenius Norm
  T FrobeniusNorm() const {
    T sum = 0;
    for (GONS_UINT i = 0; i < R; ++i) {
      for (GONS_UINT j = 0; j < C; ++j) {
        sum += (*this)(i, j) * (*this)(i, j);
      }
    }
    return std::sqrt(sum);
  }
  // 2-norm
  T Norm2() const {
    T sum = 0;
    for (GONS_UINT i = 0; i < R; ++i) {
      for (GONS_UINT j = 0; j < C; ++j) {
        sum += (*this)(i, j) * (*this)(i, j);
      }
    }
    return std::sqrt(sum);
  }

#if (ENABLE_PRINT_MATRIX == 1)
public:
  void Print() const {
    for (GONS_UINT i = 0; i < R; ++i) {
      for (GONS_UINT j = 0; j < C; ++j) {
        SHW((*this)(i, j)) << " ";
      }
      std::cout << "\n";
    }
  }
  void Print(std::string Prefix) const {
    SHW(Prefix);
    Print();
  }
#endif
#if (EMBEDDED_MODE == 0)
  friend std::ostream &operator<<(std::ostream &os,
                                  const Matrix<T, R, C> &mat) {
    for (std::size_t i = 0; i < mat.rows(); ++i) {
      for (std::size_t j = 0; j < mat.cols(); ++j) {
        os << mat(i, j) << ' ';
      }
      os << '\n';
    }
    return os;
  }
#endif


  // matrix multiply vector
  // need to be done 
  Vector<T, R> operator*(const Vector<T, C>&vec) const{
    CHECK(vec.isRowBased(), "Vector must be column based");
    Vector<T, R> result;
    for (GONS_UINT i = 0; i < R; i++) {
      for (GONS_UINT j = 0; j < C; j++) {
        result(i) += (*this)(i, j) * vec(j);
      }
    }
    return result;
  }
}; // class Matrix End

} // namespace gons
#endif