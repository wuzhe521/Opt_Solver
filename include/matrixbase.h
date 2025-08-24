/*
 *
 *
 */
#ifndef GONS_MATRIXBASE_H
#define GONS_MATRIXBASE_H

#include <iostream>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

#include <exception>

namespace gons {

template <typename T> class matrixbase {
  using string = std::string;

private:
  T **data_;
  size_t row_, col_;

private:
  void initialize(const T &num) {
    for (size_t i = 0; i < row_; ++i) {
      for (size_t j = 0; j < col_; ++j) {
        data_[i][j] = num;
      }
    }
  }

public:
  explicit matrixbase() {
    data_ = nullptr;
    row_ = col_ = 0;
  };
  matrixbase(size_t row, size_t col, const T &initValue = static_cast<T>(0))
      : row_(row), col_(col) {
    data_ = new T *[row];
    data_[0] = new T[row * col];
    for (size_t i = 1; i < row; ++i) {
      data_[i] = data_[i - 1] + col;
    }
    initialize(initValue);
  }
  matrixbase(const matrixbase &mat) {
    if (&mat == this)
      return;
    if (mat.row_ != this.row_ || mat.col_ != this.col_)
      throw std::out_of_range("diamention exceed!");
    for (size_t i = 0; i < this.row_; ++i) {
      for (size_t j = 0; j < col_; ++j) {
        data_[i][j] = mat.data_p[i][j];
      }
    }
  }

  ~matrixbase() {
    delete[] data_[0];
    delete[] data_;
  }

  std::size_t row() const { return row_; }
  std::size_t col() const { return col_; }

public:
  T &operator()(std::size_t row, std::size_t col) {
    if (row > this->row_ || col > this->col_)
      throw std::out_of_range("out of range!");
    return data_[row][col];
  }
  matrixbase &operator=(const matrixbase &mat) {
    if (this == &mat) {
      return *this;
    }

    this->~matrixbase();
    row_ = mat.row_;
    col_ = mat.col_;
    for (size_t i = 0; i < mat.row_; i++) {
      for (size_t j = 0; j < mat.col_; j++) {
        data_[i][j] = mat.data_[i][j];
      }
    }
  }
  matrixbase &operator+=(const matrixbase &mat) {
    if (mat.row_ != this->row_ || mat.col_ != this->col_)
      throw std::runtime_error("mismactch!");
    for (size_t i = 0; i < this->row_; ++i) {
      for (size_t j = 0; j < this->col_; ++j)
        data_[i][j] += mat.data_[i][j];
    }
    return *this;
  }
  matrixbase &operator+=(const T &scalar) {
    if (this->col_ == 0 || this->row_ == 0 || this->data_ == nullptr)
      throw std::runtime_error("empty matrix!");
    for (size_t i = 0; i < this->row_; ++i) {
      for (size_t j = 0; j < this->col_; ++j)
        data_[i][j] += scalar;
    }
    return *this;
  }

  matrixbase &operator-=(const matrixbase &mat) {
    if (mat.row_ != this->row_ || mat.col_ != this->col_)
      throw std::runtime_error("mismactch!");
    for (size_t i = 0; i < this->row_; ++i) {
      for (size_t j = 0; j < this->col_; ++j)
        data_[i][j] -= mat.data_[i][j];
    }
    return *this;
  }

  matrixbase &operator-=(const T &scalar) {
    if (this->col_ == 0 || this->row_ == 0 || this->data_ == nullptr)
      throw std::runtime_error("empty matrix!");
    for (size_t i = 0; i < this->row_; ++i) {
      for (size_t j = 0; j < this->col_; ++j)
        data_[i][j] -= scalar;
    }
    return *this;
  }

  matrixbase &operator*=(const matrixbase &mat) {
    if (this->col_ != mat.row_)
      throw std::runtime_error(" diamention mismatch!");
    matrixbase temp(row_, mat.cols_);
    for (int i = 0; i < temp.row_; ++i) {
      for (int j = 0; j < temp.col_; ++j) {
        for (int k = 0; k < col_; ++k) {
          temp.data_[i][j] += (data_[i][k] * mat.p[k][j]);
        }
      }
      return (*this = temp); // call overload operator= function
    }
  }
  matrixbase &operator*=(const T &scalar) {
    if (this->col_ == 0 || this->row_ == 0 || this->data_ == nullptr)
      throw std::runtime_error("empty matrix!");
    for (size_t i = 0; i < row_; ++i) {
      for (size_t j = 0; j < col_; ++j) {
        data_[i][j] *= scalar;
      }
    }
    return *this;
  }

  //matrixbase &operator/=(const matrixbase &mat) {}
  matrixbase &operator/=(const T &scalar) {
    if (this->col_ == 0 || this->row_ == 0 || this->data_ == nullptr)
      throw std::runtime_error("empty matrix!");
    for (size_t i = 0; i < row_; ++i) {
      for (size_t j = 0; j < col_; ++j) {
        data_[i][j] /= scalar;
      }
    }
    return *this;
  }
  /*
  matrixbase &operator^=(const T &scalar);
*/
public:
  friend std::ostream &operator<<(std::ostream &os, const matrixbase &mat) {
    for (size_t i = 0; i < mat.row_; ++i) {
      for (size_t j = 0; j < mat.col_; ++j) {
        os << mat.data_[i][j] << ' ';
      }
      os << '\n';
    }
    return os;
  }
};

} // namespace gons

#endif