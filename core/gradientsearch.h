#ifndef GONS_GRADIENTSEARCH_H_
#define GONS_GRADIENTSEARCH_H_

#include "function.h"
#include "linearsearch.h"
#include "matrix.h"
#include "utilites.h"
#include "vector.h"

namespace gons {
namespace gradientsearch {

using namespace gons::utilites::LOG_MSG;

// 梯度下降搜索
// 梯度采用函数的偏导数
// 步长采用 线搜索 Armijo 搜索
template <typename Function, typename X> class GradientDescentSearch {
  struct gradient_descent_parameters {
    double tolerance = 1.0e-6;
    GONS_UINT max_iterations = 1000;
    bool verbose = false;
  };
  enum class gradient_descent_state { SUCCESS, FAILURE, MAX_ITERATION_REACHED };

public:
  // Corrected constructor with matching initialization order
  GradientDescentSearch(const Function &f, const X &x)
      : f_(f), x_(x), linear_search_(f, x) {}

  ~GradientDescentSearch() = default;
  void set_params(const gradient_descent_parameters &parameters) {
    parameters_ = parameters;
  }
  void set_linear_search_params(
      const typename linearsearch::ArmijoSearch<Function, X>::ArmijoParameters &parameters) {
    linear_search_.set_params(parameters);
  }
  double SearchStep(const X &x, const Function &f) {
    return linear_search_.Search(x, f);
  }
  gradient_descent_state Optimize() {
    GONS_UINT iter = 0;
    GONS_UINT max_iter = parameters_.max_iterations;
    GONS_FLOAT tolerance = parameters_.tolerance;

    while (iter < max_iter) {
      // 计算梯度
      X gradient = f_.gradient(x_);

      // 计算步长
      double step = SearchStep(x_, f_);

      // 更新x
      X x_new = x_ - step * gradient;
      if (parameters_.verbose) // 打印信息
      {
        LOG("Iteration: " << iter);
        LOG("x = " << x_);
        LOG("f(x) = " << f_(x_));
        LOG("Step = " << step);
      }
      if (std::abs(f_(x_new) - f_(x_)) < tolerance) {
        LOG_WARNING("After " << iter << " iterations");
        LOG_WARNING("Gradient descent converged.");
        LOG_WARNING("x = " << x_);
        LOG_WARNING("f(x) = " << f_(x_));
        return gradient_descent_state::SUCCESS;
      }

      x_ = x_new;
      iter++;
    }
    LOG_WARNING("After " << iter << " iterations");
    LOG_WARNING("Gradient descent did not converge.");
    LOG_WARNING("x = " << x_);
    LOG_WARNING("f(x) = " << f_(x_));
    return gradient_descent_state::MAX_ITERATION_REACHED;
  }

private:
  Function f_;
  X x_;
  linearsearch::ArmijoSearch<Function, X> linear_search_;
  gradient_descent_parameters parameters_;
};

// Barzilai-Borwein 方法
// 参考: https://en.wikipedia.org/wiki/Barzilai%E2%80%93Borwein_method
// 参考:
// http://faculty.bicmr.pku.edu.cn/~wenzw/courses/WenyuSun_YaxiangYuan_BB.pdf

template <typename Function, typename X> class BarzilaiBorwein {
  enum class SearchMethod { BB1, BB2 };
  struct BarzilaiBorweinParameters {
    // linear search rule parameters
    double c = 0.1;
    double beta = 0.333;
    // limit step size in a reasonable range
    double alpha_upper = 1e10;
    double alpha_lower = 1e-10;
    // initial step size
    double alpha = 0.01;
    // gradient descent parameters
    double tolerance = 1.0e-6;
    int max_iterations = 10;

    bool verbose = false;
    SearchMethod method = SearchMethod::BB2;
  };
  enum class BarzilaiBorweinStatus { SUCCESS, FAILURE, MAX_ITERATION_REACHED };

public:
  BarzilaiBorwein(const Function &f, const X &x)
      : f_(f), lastGradient_(x), lastX_(x) {
    lastGradient_ = f_.gradient(lastX_);
  }
  double SearchStep(const X &x_new, const X &grad_new, const Function &f) {

    X s_k1 = x_new - lastX_;
    X y_k1 = grad_new - lastGradient_;

    double alpha = 1.0;

    switch (parameters_.method) {
    case SearchMethod::BB1: {
      // 检查分母是否为零
      FLT_EQUAL_ZERO(y_k1.dot(y_k1));
      alpha = s_k1.dot(y_k1) / y_k1.dot(y_k1);
      break;
    }
    case SearchMethod::BB2: {
      // 检查分母是否为零
      FLT_EQUAL_ZERO(y_k1.dot(y_k1));
      alpha = s_k1.dot(s_k1) / s_k1.dot(y_k1);
      break;
    }
    default:
      LOG_WARNING("Unsupported search method.");
      return 0.0;
    }
    return alpha;
  }
  X get_x() const{
    return lastX_;
  }
  // 优化方法
  BarzilaiBorweinStatus Optimize() {
    int iter = 0;

    while (iter < parameters_.max_iterations) {
      //
      X x_new = lastX_ - parameters_.alpha * lastGradient_;
      // 计算梯度
      X gradient_new = f_.gradient(x_new);

      // 更新x
      double step = SearchStep(x_new, gradient_new, f_);

      // 限制步长在合理范围内
      parameters_.alpha =
          step > parameters_.alpha_upper ? parameters_.alpha_upper : step;
      parameters_.alpha =
          step < parameters_.alpha_lower ? parameters_.alpha_lower : step;

      // 更新x 和梯度
      lastX_ = x_new;
      lastGradient_ = gradient_new;

      iter++;
      if (parameters_.verbose) // 打印信息
      {
        LOG("Iteration: " << iter);
        LOG("x = " << lastX_);
        LOG("f(x) = " << f_(lastX_));
        LOG("Step = " << step);
      }

      // 检查梯度是否足够小
      if (lastGradient_.Norm2() < parameters_.tolerance) {
        LOG_WARNING("After " << iter << " iterations");
        LOG_WARNING("Barzilai-Borwein converged.");
        LOG_WARNING("x = " << lastX_);
        LOG_WARNING("f(x) = " << f_(lastX_));
        return BarzilaiBorweinStatus::SUCCESS;
      }
    }

    // 检查是否达到最大迭代次数
    if (iter == parameters_.max_iterations) {
      LOG_WARNING("After " << iter << " iterations");
      LOG_WARNING("Barzilai-Borwein did not converge.");
      LOG_WARNING("x = " << lastX_);
      LOG_WARNING("f(x) = " << f_(lastX_));
      return BarzilaiBorweinStatus::MAX_ITERATION_REACHED;
    }
  }

private:
  Function f_;
  X lastGradient_;
  X lastX_;
  BarzilaiBorweinParameters parameters_;
};
} // namespace gradientsearch
} // namespace gons

#endif // GONS_GRADIENTSEARCH_H_