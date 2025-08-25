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

  matrixbase(const matrixbase<T> &mat) {
    if (&mat == this)
      return;
    if (data_ != nullptr)
      this->~matrixbase();
    this->row_ = mat.row_;
    this->col_ = mat.col_;
    data_ = new T *[mat.row_];
    data_[0] = new T[mat.row_ * mat.col_];
    for (size_t i = 1; i < mat.row_; ++i) {
      data_[i] = data_[i - 1] + mat.col_;
    }
    for (size_t i = 0; i < mat.row_; ++i) {
      for (size_t j = 0; j < mat.col_; ++j) {
        this->data_[i][j] = mat.data_[i][j];
      }
    }
  }
  ~matrixbase() {
    delete[] data_[0];
    delete[] data_;
  }

  std::size_t rows() const { return row_; }
  std::size_t cols() const { return col_; }

public:
  T &operator()(std::size_t row, std::size_t col) {
    if (row > this->row_ || col > this->col_)
      throw std::out_of_range("out of range!");
    return data_[row][col];
  }
  matrixbase<T> &operator=(const matrixbase<T> &mat) {
    if (this == &mat) {
      return *this;
    }
    this->~matrixbase();
    row_ = mat.row_;
    col_ = mat.col_;
    data_ = new T *[mat.row_];
    data_[0] = new T[mat.row_ * mat.col_];
    for (size_t i = 1; i < mat.row_; ++i) {
      data_[i] = data_[i - 1] + mat.col_;
    }
    for (size_t i = 0; i < mat.row_; i++) {
      for (size_t j = 0; j < mat.col_; j++) {
        data_[i][j] = mat.data_[i][j];
      }
    }
  }
  matrixbase<T> &operator+=(const matrixbase<T> &mat) {
    if (mat.row_ != this->row_ || mat.col_ != this->col_)
      throw std::runtime_error("mismactch!");
    for (size_t i = 0; i < this->row_; ++i) {
      for (size_t j = 0; j < this->col_; ++j)
        data_[i][j] += mat.data_[i][j];
    }
    return *this;
  }
  matrixbase<T> &operator+=(const T &scalar) {
    if (this->col_ == 0 || this->row_ == 0 || this->data_ == nullptr)
      throw std::runtime_error("empty matrix!");
    for (size_t i = 0; i < this->row_; ++i) {
      for (size_t j = 0; j < this->col_; ++j)
        data_[i][j] += scalar;
    }
    return *this;
  }

  matrixbase<T> &operator-=(const matrixbase<T> &mat) {
    if (mat.row_ != this->row_ || mat.col_ != this->col_)
      throw std::runtime_error("mismactch!");
    for (size_t i = 0; i < this->row_; ++i) {
      for (size_t j = 0; j < this->col_; ++j)
        data_[i][j] -= mat.data_[i][j];
    }
    return *this;
  }

  matrixbase<T> &operator-=(const T &scalar) {
    if (this->col_ == 0 || this->row_ == 0 || this->data_ == nullptr)
      throw std::runtime_error("empty matrix!");
    for (size_t i = 0; i < this->row_; ++i) {
      for (size_t j = 0; j < this->col_; ++j)
        data_[i][j] -= scalar;
    }
    return *this;
  }

  matrixbase<T> &operator*=(const matrixbase<T> &mat) {
    if (this->col_ != mat.row_)
      throw std::runtime_error("dimension mismatch!");
    matrixbase<T> temp(row_, mat.cols_);
    for (int i = 0; i < temp.row_; ++i) {
      for (int j = 0; j < temp.col_; ++j) {
        for (int k = 0; k < col_; ++k) {
          temp.data_[i][j] += (data_[i][k] * mat.data_[k][j]);
        }
      }
      return (*this = temp); // call overload operator= function
    }
  }
  matrixbase<T> &operator*=(const T &scalar) {
    if (this->col_ == 0 || this->row_ == 0 || this->data_ == nullptr)
      throw std::runtime_error("empty matrix!");
    for (size_t i = 0; i < row_; ++i) {
      for (size_t j = 0; j < col_; ++j) {
        data_[i][j] *= scalar;
      }
    }
    return *this;
  }

  // matrixbase &operator/=(const matrixbase &mat) {}
  matrixbase<T> &operator/=(const T &scalar) {
    if (this->col_ == 0 || this->row_ == 0 || this->data_ == nullptr)
      throw std::runtime_error("empty matrix!");
    for (std::size_t i = 0; i < row_; ++i) {
      for (std::size_t j = 0; j < col_; ++j) {
        data_[i][j] /= scalar;
      }
    }
    return *this;
  }
  static matrixbase<T> Indetity(std::size_t size) {
    matrixbase<T> idet(size, size);
    for (std::size_t i = 0; i < size; ++i) {
      idet.data_[i][i] = T(1);
    }
    return idet;
  }
  matrixbase<T> inverse() {
    if (this->row_ != this->col_)
      throw std::runtime_error("not square matrix!");
    matrixbase<T> matI = Indetity(this->row_);
    matrixbase<T> oriI(*this);
    auto vectorScale = [this](T *p, std::size_t len, T &scalar) {
      for (std::size_t idx = 0; idx < len; ++idx) {
        p[idx] *= scalar;
      }
    };
    auto vectorSwap = [this](T *p, T *l) {
      for (std::size_t idx = 0; idx < len; ++idx) {
        T temp = p[idx], p[idx] = l[idx], l[idx] = temp;
      }
    };
    auto vectoradd = [this](T *p, T *l) {
      for (std::size_t idx = 0; idx < len; ++idx) {
        p[idx] += l[idx];
      }
    };

    for (std::size_t i = 0; i < this->col_; ++i) {
      std::size_t loc = 0;
      for (std::size_t j = 1; j < this->row_; ++j) {
        if (data_[j][i] > data_[loc][i]) {
          loc = j;
        }
      }
      vectorSwap(data_[0], data_[loc]);
      vectorSwap(matI.data_[0], matI.data_[loc]);
    }

    return oriI;
  }

  matrixbase &operator^=(const T &scalar);

public:
  friend std::ostream &operator<<(std::ostream &os, const matrixbase<T> &mat) {
    for (std::size_t i = 0; i < mat.row_; ++i) {
      for (std::size_t j = 0; j < mat.col_; ++j) {
        os << mat.data_[i][j] << ' ';
      }
      os << '\n';
    }
    return os;
  }
  friend matrixbase<T> operator+(const matrixbase<T> &matA,
                                 const matrixbase<T> &matB) {
    if (matA.row_ != matB.row_ || matA.col_ != matB.col_)
      throw std::runtime_error("dimension unequal!");
    matrixbase<T> result(matA.row_, matA.col_);
    for (std::size_t i = 0; i < matA.row_; ++i) {
      for (std::size_t j = 0; j < matA.col_; ++j) {
        result.data_[i][j] = matA.data_[i][j] + matB.data_[i][j];
      }
    }
    return result;
  }
  friend matrixbase<T> operator-(const matrixbase<T> &matA,
                                 const matrixbase<T> &matB) {
    if (matA.row_ != matB.row_ || matA.col_ != matB.col_)
      throw std::runtime_error("dimension unequal!");
    matrixbase<T> result(matA.row_, matA.col_);
    for (std::size_t i = 0; i < matA.row_; ++i) {
      for (std::size_t j = 0; j < matA.col_; ++j) {
        result.data_[i][j] = matA.data_[i][j] - matB.data_[i][j];
      }
    }
    return result;
  }
  friend matrixbase<T> operator*(const matrixbase<T> &matA,
                                 const matrixbase<T> &matB) {
    if (matA.col_ != matB.row_)
      throw std::runtime_error("dimension mismatch!");
    matrixbase<T> result(matA.row_, matB.col_);
    for (std::size_t i = 0; i < result.row_; ++i) {
      for (std::size_t j = 0; j < result.col_; ++j) {
        for (std::size_t k = 0; k < matA.col_; ++k) {
          result.data_[i][j] += (matA.data_[i][k] * matB.data_[k][j]);
        }
      }
    }
    return result;
  }
};

} // namespace gons

#endif
