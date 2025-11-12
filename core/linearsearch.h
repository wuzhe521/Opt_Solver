#ifndef GONS_LINEARSEARCH_H
#define GONS_LINEARSEARCH_H
#include "function.h"
#include "matrix.h"
#include "utilites.h"
#include "vector.h"

namespace gons {
namespace linearsearch {
using namespace gons::utilites::LOG_MSG;
template <typename Function, typename X> class ArmijoSearch {

public:
  struct ArmijoParameters {
    double alpha = 1e-4;
    double beta = 0.2;
    double gamma = 0.333;
    double epsilon = GONS_FLT_EPSILON;
    GONS_UINT max_inner_iter = 1000u;
    bool enable_max_iter = true;
    GONS_UINT max_iter = 1000u;
    bool print_info = false;
  };

public:
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
      X x_new = x_ - beta * gradient;
      double f_x_new = f_(x_new);
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
  // in order to use the same search rule for different function, we need to
  // pass the function and initial point as arguments
  double Search(const X &x0, const Function &f) {
    const double min_beta =
        1e-10; // Minimum alpha to avoid numerical instability
    GONS_UINT iteration = 0;
    double beta = params_.beta;
    while (iteration < params_.max_inner_iter) {

      // Check if alpha is too small
      if (params_.beta < min_beta) {
        break;
      }

      // Cache gradient to avoid repeated computation
      X gradient = f.gradient(x0);
      double f_x = f(x0);
      double f_x_new = f(x0 - beta * gradient);
      if (f_x_new < f_x - params_.alpha * beta * gradient.Norm2()) {
        // Armijo condition is satisfied
        params_.beta = beta;
        return beta;
      }
      // Update alpha for next iteration
      beta *= params_.gamma;
      if (params_.enable_max_iter)
        iteration++;
    }
    params_.beta = beta;
    return beta;
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

// Goldstein search
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

public:
  enum class GoldsteinStatus { Success, Failure };

public:
  GoldsteinSearch(const Function &f, const X &x0) : f_(f), x_(x0) {}
  void set_params(const GoldsteinParameters &params) { params_ = params; }
  X get_x() const { return x_; }
  double get_function_value() const { return f_(x_); }

  X GoldsteinRule() {
    // Implementation of Goldstein search
    const double min_beta =
        1e-10;                   // Minimum beta to avoid numerical instability
    const double max_beta = 1.0; // Maximum beta to avoid divergence
    GONS_UINT iteration = 0;
    double beta = params_.beta;
    X gradient = f_.gradient(x_);
    double f_x = f_(x_);
    double gradient_norm = gradient.Norm2();

    while (true) {
      ++iteration;
      if (iteration > params_.max_inner_iter || beta < min_beta ||
          beta > max_beta) {
        break;
      }

      double f_x_new = f_(x_ - beta * gradient);
      if (f_x_new < f_x - beta * params_.alpha * gradient_norm &&
          f_x_new > f_x - (1 - params_.alpha) * beta * gradient_norm) {
        params_.beta = beta;
        return x_ - params_.alpha * gradient;
      }
      beta *= params_.gamma;
    }
    return x_ - params_.alpha * gradient;
  }
  double Search(const X &x0, const Function &f) {
    // Implementation of Goldstein search
    const double min_beta =
        1e-10;                   // Minimum beta to avoid numerical instability
    const double max_beta = 1.0; // Maximum beta to avoid divergence
    GONS_UINT iteration = 0;
    double beta = params_.beta;
    X gradient = f.gradient(x0);
    double f_x = f(x0);
    double gradient_norm = gradient.Norm2();

    while (true) {
      ++iteration;
      if (iteration > params_.max_inner_iter || beta < min_beta ||
          beta > max_beta) {
        break;
      }
      double f_x_new = f(x0 - beta * gradient);
      if (f_x_new < f_x - beta * params_.alpha * gradient_norm &&
          f_x_new > f_x - (1 - params_.alpha) * beta * gradient_norm) {
        params_.beta = beta;
        return beta;
      }
      beta *= params_.gamma;
    }
    return beta;
  }
  GoldsteinStatus Optimize() {
    LOG(SOLVER_HEADER);
    GONS_UINT i = 0;
    while (true) {
      ++i;
      if (params_.enable_max_iter)
        if (i >= params_.max_iter)
          break;
      X x_new = GoldsteinRule();
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
};

// wolfe search

template <typename Function, typename X> class WolfeSearch {
public:
  struct WolfeParameters {
    double alpha_1 = 0.0001;
    double alpha_2 = 0.9;
    double beta = 0.2;
    double gamma = 0.333;
    double epsilon = GONS_FLT_EPSILON;
    GONS_UINT max_inner_iter = 1000u;
    bool enable_max_iter = true;
    GONS_UINT max_iter = 1000u;
    bool print_info = false;
  };

public:
  enum class WolfeStatus { Success, Failure };

public:
  WolfeSearch(Function f, X x) : f_(f), x_(x) {}
  void set_parameters(const WolfeParameters &params) { params_ = params; }
  X get_x() const { return x_; }
  double get_function_value() const { return f_(x_); }
  X WolfeRule() {
    const double min_beta = 1e-10;
    GONS_UINT iteration = 0;
    double beta = params_.beta;
    while (true) {
      ++iteration;
      if (iteration > params_.max_inner_iter) {
        break;
      }
      if (beta < min_beta) { // Check if beta is too small and break if it is
        break;
      }
      X gradient = f_.gradient(x_);
      double f_x = f_(x_);
      // Wolfe conditions
      double f_x_new = f_(x_ - beta * gradient);
      if (f_x_new <= f_x - params_.alpha_1 * beta * gradient.Norm2()) {
        // Armijo condition satisfied
        X gradient_new = f_.gradient(x_ - beta * gradient);
        if (-gradient_new.Norm2() > -params_.alpha_2 * gradient.Norm2()) {
          // Wolfe condition satisfied
          params_.beta = beta;
          return x_ - beta * gradient;
        } else {
          // Only Armijo condition satisfied, decrease beta
          beta *= params_.gamma;
        }
      } else {
        // Neither condition satisfied, increase beta
        beta *= params_.gamma * 1.1;
      }
    }
    params_.beta = beta;
    return x_ - beta * f_.gradient(x_);
  }

  double Search(Function &f, X &d, X &x) {
    const double min_beta = 1e-10;
    GONS_UINT iteration = 0;
    double beta = params_.beta;
    while (true) {
      ++iteration;
      if (iteration > params_.max_inner_iter) {
        break;
      }
      if (beta < min_beta) { // Check if beta is too small and break if it is
        break;
      }
      X gradient = f.gradient(x);
      double f_x = f(x);
      // Wolfe conditions
      double f_x_new = f(x - beta * gradient);
      if (f_x_new <= f_x - params_.alpha_1 * beta * gradient.Norm2()) {
        // Armijo condition satisfied
        X gradient_new = f.gradient(x - beta * gradient);
        if (-gradient_new.Norm2() > -params_.alpha_2 * gradient.Norm2()) {
          // Wolfe condition satisfied
          params_.beta = beta;
          return beta;
        } else {
          // Only Armijo condition satisfied, decrease beta
          beta *= params_.gamma;
        }
      } else {
        // Armijo condition not satisfied, increase beta
        beta *= params_.gamma * 1.1;
      }
    }
    return beta;
  }
  WolfeStatus Optimize() {
    LOG(SOLVER_HEADER);
    GONS_UINT i = 0;
    while (true) {
      ++i;
      if (params_.enable_max_iter)
        if (i >= params_.max_iter)
          break;
      X x_new = WolfeRule();
      if (params_.print_info) {
        LOG("Iteration: " << i);
        LOG_WARNING("Function value last: " << f_(x_));
        LOG_WARNING("Function value new: " << f_(x_new));
        LOG_WARNING("Gradient Norm2: " << f_.gradient(x_).Norm2());
        LOG_ERROR("x: " << x_);
      }
      if (std::abs(f_(x_new) - f_(x_)) < params_.epsilon) {
        LOG_ERROR("After " << i << " iterations, ");
        LOG_ERROR("Wolfe search converged to a local minimum");
        return WolfeStatus::Success;
      }
      x_ = x_new;
    }
    return WolfeStatus::Failure;
  }

private:
  Function f_;
  X x_;
  //
  WolfeParameters params_;
};
} // namespace linearsearch
} //   namespace gons
#endif