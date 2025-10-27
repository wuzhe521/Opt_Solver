#ifndef GONS_LINEARSEARCH_H
#define GONS_LINEARSEARCH_H
#include "function.h"
#include "matrix.h"
#include "utilites.h"
#include "vector.h"

namespace gons {

using namespace gons::utilites::LOG_MSG;

template <typename Function, typename X> class ArmijoSearch {
public:
  struct ArmijoParameters {
    double alpha = 1e-4;
    double beta = 0.5;
    double gamma = 0.333;
    double epsilon = GONS_FLT_EPSILON;
    GONS_UINT max_inner_iter = 1000u;
    bool enable_max_iter = true;
    GONS_UINT max_iter = 1000u;
    bool print_info = false;
  };

  enum class ArmijoStatus { Success, Failure };

public:
  ArmijoSearch(const Function &f, const X &x) : f_(f), x_(x) {}
  void set_params(const ArmijoParameters &params) { params_ = params; }
  X get_x() const { return x_; }
  double get_function_value() const { return f_(x_); }
  X ArmijoRule() {
    const double min_beta =
        1e-10; // Minimum alpha to avoid numerical instability

    GONS_UINT iteration = 0;

    while (iteration < params_.max_inner_iter) {
      double beta = params_.beta;
      // Check if alpha is too small
      if (params_.beta < min_beta) {
        break;
      }

      // Cache gradient to avoid repeated computation
      X gradient = f_.gradient(x_);
      double f_x = f_(x_);
      double f_x_new = f_(x_ - beta * gradient);
      if (f_x_new < f_x - params_.alpha * beta * gradient.Norm2()) {
        // Armijo condition is satisfied
        params_.beta = beta;
        return x_ - beta * gradient;
      }
      // Update alpha for next iteration
      beta *= params_.gamma;
      if (params_.enable_max_iter)
        iteration++;
    }
    return x_ - params_.beta * f_.gradient(x_);
  }
  ArmijoStatus Optimize() {
    LOG(SOLVER_HEADER);
    GONS_UINT i = 0;
    while (true) {
      ++i;
      if (params_.enable_max_iter && i >= params_.max_iter) {
        LOG_ERROR("Max iterations reached: " << params_.max_iter);
        break;
      }
      X x_new = ArmijoRule();
      double f_old = f_(x_);
      double f_new = f_(x_new);
      double grad_norm = f_.gradient(x_).Norm2();
      if (params_.print_info) {
        LOG("Iteration: " << i);
        LOG_WARNING("Function value last: " << f_old);
        LOG_WARNING("Function value new: " << f_new);
        LOG_WARNING("Gradient Norm2: " << grad_norm);
        LOG_ERROR("x: " << x_);
      }
      if (std::abs(f_new - f_old) < params_.epsilon) {
        LOG_ERROR("After " << i << " iterations, ");
        LOG_ERROR("Armijo search converged to a local minimum");
        return ArmijoStatus::Success;
      }
      x_ = x_new;
    }
    return ArmijoStatus::Failure;
  }
  //
private:
  Function f_;
  X x_;

  //  Paremeters
  ArmijoParameters params_;
};
/*
template <typename Function, typename X> class GoldsteinSearch {
public:
  struct GoldsteinParameters {
    double alpha = 0.25;
    double beta = 0.5;
    double gamma = 0.333;
    double epsilon = GONS_FLT_EPSILON;
    GONS_UINT max_inner_iter = 1000u;
    bool enable_max_iter = true;
    GONS_UINT max_iter = 1000u;
    bool print_info = false;
  };
  enum class GoldsteinStatus { Success, Failure };

public:
  GoldsteinSearch(const Function &f, const X &x0,
                  const GoldsteinParameters &params)
      : f_(f), x_(x0), params_(params) {}

  X get_x() const { return x_; }
  double get_function_value() const { return f_(x_); }

  X Search() {
    // Implementation of Goldstein search
    const double min_alpha =
        1e-10; // Minimum alpha to avoid numerical instability]
    GONS_UINT iteration = 0;
    double alpha = params_.alpha;
    while (true) {
      //
      X gradient = f_.gradient(x_);
      double f_x = f_(x_);
      double f_x_new = f_(x_ - params_.alpha * gradient);
      if (f_x_new <=
          f_x - (1 - params_.alpha) * params_.beta * gradient.Norm2()) {
        //
        alpha *= (1 + params_.gamma);
      }
      if (f_x_new >= f_x - params_.alpha * params_.beta * gradient.Norm2()) {
        alpha *= params_.gamma;
      }
      if (f_x_new < f_x - params_.beta * alpha * gradient.Norm2() &&
          f_x_new > f_x - (1 - alpha) * params_.beta * gradient.Norm2()) {
        break;
      }
      if (alpha < min_alpha) {
        break;
      }
      x_ -= params_.alpha * gradient;
    }
  }
  GoldsteinStatus Optimize() {
    LOG(SOLVER_HEADER);
    GONS_UINT i = 0;
    while (true) {
      ++i;
      if (params_.enable_max_iter)
        if (i >= params_.max_iter)
          break;
      X x_new = Search();
      if (params_.print_info) {
        LOG("Iteration: " << i);
        LOG_WARNING("Function value last: " << f_(x_));
        LOG_WARNING("Function value new: " << f_(x_new));
        LOG_WARNING("Gradient Norm2: " << f_.gradient(x_).Norm2());
        LOG_ERROR("x: " << x_);
      }
      if (std::abs(f_(x_new) - f_(x_)) < params_.epsilon) {
        LOG_ERROR("After " << i << " iterations, ");
        LOG_ERROR("Goldstein search converged to a local minimum");
        return GoldsteinStatus::Success;
      }
      x_ = x_new;
    }
    return GoldsteinStatus::Failure;
  }

private:
  Function f_;
  X x_;

  // Parameters
  GoldsteinParameters params_;
};*/

} //   namespace gons
#endif