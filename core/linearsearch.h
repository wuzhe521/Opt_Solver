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
    const double min_alpha =
        1e-10; // Minimum alpha to avoid numerical instability

    int iteration = 0;

    while (iteration < params_.max_inner_iter) {
      // Check if alpha is too small
      if (params_.alpha < min_alpha) {
        break;
      }

      // Cache gradient to avoid repeated computation
      X gradient = f_.gradient(x_);
      double f_x = f_(x_);
      double f_x_new = f_(x_ - params_.alpha * gradient);
      if (f_x_new < f_x - params_.beta * params_.alpha * gradient.Norm2()) {
        // Armijo condition is satisfied
        return x_ - params_.alpha * gradient;
      }
      // Update alpha for next iteration
      params_.alpha *= params_.gamma;
      if (params_.enable_max_iter)
        iteration++;
    }
    return x_ - params_.alpha * f_.gradient(x_);
  }
  ArmijoStatus Optimize() {
    LOG(SOLVER_HEADER);
    GONS_UINT i = 0;
    while(true) {
      ++i;
      if (params_.enable_max_iter)
          if (i >= params_.max_iter) break;
      X x_new = ArmijoRule();
      if (params_.print_info) {
        LOG("Iteration: " << i);
        LOG_WARNING("Function value last: " << f_(x_));
        LOG_WARNING("Function value new: " << f_(x_new));
        LOG_WARNING("Gradient Norm2: " << f_.gradient(x_).Norm2());
        LOG_ERROR("x: " << x_);
      }
      if (std::abs(f_(x_new) - f_(x_)) < params_.epsilon) {
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

} // namespace gons
#endif //