#ifndef GONS_CORE_VECTOR_H_
#define GONS_CORE_VECTOR_H_

#include "config.h"
#include "matrix.h"
#include "utilites.h"

namespace gons {

template <typename T, GONS_UINT R, GONS_UINT C> class Matrix;

template <typename T, GONS_UINT N> class Vector : public Matrix<T, 1, N> {
  using Matrix<T, 1, N>::Matrix;
protected:
  bool row_based = true;
public:
  bool getRowBased() const { return this->row_based; }
  void setRowBased(bool row_based) { this->row_based = row_based; }
public:
  Vector() = default;
  Vector(const Matrix<T, 1, N> &m) : Matrix<T, 1, N>(m) {}
  Vector(const std::initializer_list<T> &list) : Matrix<T, 1, N>(list) {

    CHECK(list.size() != N, "list size not equal to N")

    GONS_UINT i = 0;
    for (auto &item : list) {
      this->data_[i][0] = item;
      i++;
      if (i >= N)
        break;
    }
  }

  Vector(const Vector<T, N> &v) {
    for (GONS_UINT i = 0; i < N; i++) {
      this->data_[0][i] = v.data_[0][i];
    }
  }

  bool isRowBased() const { return row_based; }
  bool isColumnBased() const { return !row_based; }
  Vector<T, N> &operator=(const Vector<T, N> &v) {
    for (GONS_UINT i = 0; i < N; i++) {
      this->data_[0][i] = v.data_[0][i];
    }
    return *this;
  }

  Vector<T, N> &operator=(const Matrix<T, N, 1> &m) {
    for (GONS_UINT i = 0; i < N; i++) {
      this->data_[0][i] = m(i, 0);
    }
    return *this;
  }
  T &operator()(GONS_UINT i) { return this->data_[0][i]; }
  const T &operator()(GONS_UINT i) const { return this->data_[0][i]; }

  T dot(const Vector<T, N> &v) const {
    T result = 0;
    for (GONS_UINT i = 0; i < N; i++) {
      result += this->data_[0][i] * v(i);
    }
    return result;
  }
  Vector<T, N> operator+(const Vector<T, N> &v) const {
    Vector<T, N> result;
    for (GONS_UINT i = 0; i < N; i++) {
      result(i) = this->data_[0][i] + v(i);
    }
    return result;
  }
  Vector<T, N> operator-(const Vector<T, N> &v) const {
    Vector<T, N> result;
    for (GONS_UINT i = 0; i < N; i++) {
      result(i) = this->data_[0][i] - v(i);
    }
    return result;
  }
  Vector<T, N> operator+=(const Vector<T, N> &v) {
    for (GONS_UINT i = 0; i < N; i++) {
      this->data_[0][i] += v(i);
    }
    return *this;
  }
  Vector<T, N> operator-=(const Vector<T, N> &v) {
    for (GONS_UINT i = 0; i < N; i++) {
      this->data_[0][i] -= v(i);
    }
    return *this;
  }
  Vector<T, N> operator*(const T &t) const {
    Vector<T, N> result;
    for (GONS_UINT i = 0; i < N; i++) {
      result(i) = this->data_[0][i] * t;
    }
    return result;
  }
  Vector<T, N> operator/(const T &t) const {
    Vector<T, N> result;
    for (GONS_UINT i = 0; i < N; i++) {
      result(i) = this->data_[0][i] / t;
    }
    return result;
  }

  // vector multiply matrix
  Vector<T,N> transpose() const {
    Vector<T,N> result;
    for (GONS_UINT i = 0; i < N; i++) {
      result(i) = this->data_[0][i];
    }
    result.setRowBased(!result.getRowBased());
    return result;
  }

  Matrix<T, N, N> outerProduct(const Vector<T, N> &v) const {
    Matrix<T, N, N> result;
    for (GONS_UINT i = 0; i < N; i++) {
      for (GONS_UINT j = 0; j < N; j++) {
        result(i, j) = this->data_[0][i] * v(j);
      }
    }
    return result;
  }
  T Sum() const {
    T result = 0;
    for (GONS_UINT i = 0; i < N; i++) {
      result += this->data_[0][i];
    }
    return result;
  }

 
};
} // namespace gons
#endif // GONS_CORE_VECTOR_H_