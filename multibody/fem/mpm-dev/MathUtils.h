#pragma once

#include "drake/common/eigen_types.h"

namespace drake {
namespace multibody {
namespace mpm {

// A class with some linear/tensor algebra helper routines
class MathUtils {
 public:
    // Return (i, j, k)th entry of the third order permutation tensor
    static double LeviCivita(int i, int j, int k);

    // Calculate A:ε
    static Vector3<double> APermutation(const Matrix3<double>& A);
};  // class MathUtils

}  // namespace mpm
}  // namespace multibody
}  // namespace drake
