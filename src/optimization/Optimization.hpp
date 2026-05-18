// Future: shape/parametric optimization driver. Interface only.
#pragma once
#include <functional>
#include <vector>
namespace simall::optimization {
struct DesignVariable { double value, lower, upper; };
using ObjectiveFn = std::function<double(const std::vector<double>&)>;
class IOptimizer {
public:
    virtual ~IOptimizer() = default;
    virtual std::vector<double> minimize(ObjectiveFn,
        std::vector<DesignVariable>, int maxIter) = 0;
};
}
