#ifndef Constants_h
#define Constants_h

#include <cmath>
#include <limits>
#include <string>

namespace OpenSees
{
namespace Constants
{
    inline constexpr double         pi = 3.141592653589793238462643383279502884L;
    inline constexpr double     inv_pi = 0.318309886183790671537767526745028724L;
    inline constexpr double        inf = std::numeric_limits<double>::infinity();
    inline const double            nan = std::nan("1");

    inline constexpr double      sqrt3 = 1.732050807568877293527446341505872367L;
    inline constexpr double  inv_sqrt3 = 0.577350269189625764509148780501957456L;

    inline constexpr double      sqrt2 = 1.414213562373095048801688724209698079L;
}
}

#endif // Constants_h
