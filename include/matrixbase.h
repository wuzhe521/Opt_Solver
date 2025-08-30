/*
 *
 *
 */
#ifndef GONS_MATRIXBASE_H
#define GONS_MATRIXBASE_H

#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

#include <exception>

// 1: ON 0: OFF
#define MATRIX_INV_MAJOR_ELE_SWT 1

namespace gons {

template <typename T> class matrixbase {
private:
  T **data_;
  std::size_t row_;
  std::size_t col_;

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
  matrixbase(std::initializer_list<std::initializer_list<T>> list) {

    auto diamension_valid = [&list]() {
      if (list.size() == 0)
        return false;
      for (auto itr = list.begin(); itr != list.end(); ++itr) {
        if (list.begin()->size() != itr->size())
          return false;
      }
      return true;
    };

    if (!diamension_valid())
      throw std::runtime_error("diamension mismactch!");
    this->row_ = list.size();
    this->col_ = list.begin()->size();
    this->data_ = new T *[this->row_];
    this->data_[0] = new T[this->row_ * this->col_];
    std::copy(list.begin()->begin(), list.begin()->end(), this->data_[0]);
    for (size_t i = 1; i < this->row_; ++i) {
      data_[i] = data_[i - 1] + this->col_;
      std::copy((list.begin() + i)->begin(), (list.begin() + i)->end(),
                this->data_[i]);
    }
  }
  ~matrixbase() {
    if (data_ != nullptr) {
      delete[] data_[0];
      delete[] data_;
    }
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
    matrixbase<T> temp(row_, mat.col_);
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

  matrixbase<T> Identity() {
    if (this->col_ != this->row_ || this->col_ == 0)
      throw std::runtime_error("not valid square matrix");
    const auto &size = this->col_;
    matrixbase<T> idet(size, size);
    for (std::size_t i = 0; i < size; ++i) {
      idet.data_[i][i] = T(1);
    }
    return idet;
  }

  static matrixbase<T> Augment(const matrixbase<T> &orig) {
    if (orig.row_ != orig.col_ || orig.row_ == 0 || orig.col_ == 0) {
      throw std::runtime_error("not valid square matrix");
    }
    matrixbase<T> aug(orig.row_, orig.col_ * 2);
    for (std::size_t i = 0; i < orig.row_; ++i) {
      for (std::size_t j = 0; j < orig.col_; ++j) {
        aug.data_[i][j] = orig.data_[i][j];
      }
    }
    for (std::size_t i = 0; i < orig.row_; ++i) {
      aug.data_[i][i + orig.col_] = T(1);
    }
    return aug;
  }
  matrixbase<T> Augment() {
    if (row_ != col_ || row_ == 0 || col_ == 0) {
      throw std::runtime_error("not valid square matrix");
    }
    matrixbase<T> aug(row_, col_ * 2);
    for (std::size_t i = 0; i < row_; ++i) {
      for (std::size_t j = 0; j < col_; ++j) {
        aug.data_[i][j] = data_[i][j];
      }
    }
    for (std::size_t i = 0; i < row_; ++i) {
      aug.data_[i][i + col_] = T(1);
    }
    return aug;
  }
  matrixbase<T> transpose() {
    if (this->row_ == 0 || this->col_ == 0)
      throw std::runtime_error("empty matrix!");
    matrixbase<T> trans(this->col_, this->row_);
    for (std::size_t i = 0; i < this->row_; ++i) {
      for (std::size_t j = 0; j < this->col_; ++j) {
        trans.data_[j][i] = this->data_[i][j];
      }
    }
    return trans;
  }
  matrixbase<T> inverse() {
    if (this->row_ != this->col_ || this->row_ < 1)
      throw std::runtime_error("not square matrix!");
    if (this->row_ == 1) {
      if (this->data_[0][0] == 0)
        throw std::runtime_error(" matrix uninvertible");
      matrixbase<T> inv(1, 1);
      inv.data_[0][0] = 1 / this->data_[0][0];
      return inv;
    }
    if (this->row_ == 2) {
      matrixbase<T> inv(2, 2);
      T det = this->data_[0][0] * this->data_[1][1] -
              this->data_[0][1] * this->data_[1][0];
      if (det == 0)
        throw std::runtime_error(" matrix uninvertible");
      inv.data_[0][0] = this->data_[1][1] / det;
      inv.data_[0][1] = -this->data_[0][1] / det;
      inv.data_[1][0] = -this->data_[1][0] / det;
      inv.data_[1][1] = this->data_[0][0] / det;
      return inv;
    }
    matrixbase<T> augMat(this->Augment());
    const auto &row_cnt = this->row_;
    const auto &col_cnt = this->col_;
    const auto &size = this->row_;
    const auto &aug_col = augMat.col_;
    auto vector_scal = [](T *p, const double &scale, const std::size_t &len) {
      for (std::size_t idx = 0; idx < len; ++idx) {
        p[idx] *= scale;
      }
    };
    auto vector_scal_sub = [this](T *p, T *l, const T &scale,
                                  const std::size_t &len) {
      for (std::size_t idx = 0; idx < len; ++idx) {
        auto temp = l[idx] * scale;
        p[idx] = p[idx] - temp;
      }
    };
    auto vector_switch = [this](T *p, T *l, const std::size_t &len) {
      T temp;
      for (std::size_t i = 0; i < len; ++i) {
        temp = p[i], p[i] = l[i], l[i] = temp;
      }
    };
    for (std::size_t i = 0; i < col_cnt; ++i) {

      std::size_t max_ele_loc = 0;
      T max_ele = T(0);
      for (std::size_t j = i; j < row_cnt; ++j) {
        if (std::abs(augMat.data_[j][i]) > max_ele) {
          max_ele = std::abs(augMat.data_[j][i]);
          max_ele_loc = j;
        }
      }
      vector_switch(augMat.data_[i], augMat.data_[max_ele_loc], aug_col);
      if (std::abs(augMat.data_[i][i]) < 1E-10)
        throw std::runtime_error(" matrix uninvertible");
      const double &para = 1 / augMat.data_[i][i];
      vector_scal(augMat.data_[i], para, aug_col);
      for (std::size_t j = 0; j < row_cnt; ++j) {
        if (i == j)
          continue;
        const double para = augMat.data_[j][i];
        vector_scal_sub(augMat.data_[j], augMat.data_[i], para, aug_col);
      }
    }
    matrixbase<T> inv(this->Identity());
    for (std::size_t i = 0; i < row_cnt; ++i) {
      for (std::size_t j = 0; j < col_cnt; ++j) {
        inv.data_[i][j] = augMat.data_[i][j + col_cnt];
      }
    }
    return inv;
  }

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
