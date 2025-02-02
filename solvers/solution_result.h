#pragma once

#include <ostream>
#include <string>

#include "drake/common/drake_deprecated.h"

namespace drake {
namespace solvers {
enum SolutionResult {
  kSolutionFound = 0,           ///< Found the optimal solution.
  kInvalidInput = -1,           ///< Invalid input.
  kInfeasibleConstraints = -2,  ///< The primal is infeasible.
  kUnbounded = -3,              ///< The primal is unbounded.
  kUnknownError = -4,           ///< Unknown error.
  kInfeasibleOrUnbounded =
      -5,                ///< The primal is either infeasible or unbounded.
  kIterationLimit = -6,  ///< Reaches the iteration limits.
  kDualInfeasible = -7,  ///< Dual problem is infeasible. In this case we cannot
                         /// infer the status of the primal problem.
  kInfeasible_Or_Unbounded
  DRAKE_DEPRECATED("2022-07-01", "Use kInfeasibleOrUnbounded instead.") = -5,
};

std::string to_string(SolutionResult solution_result);
std::ostream& operator<<(std::ostream& os, SolutionResult solution_result);
}  // namespace solvers
}  // namespace drake
