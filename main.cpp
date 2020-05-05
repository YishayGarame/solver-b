#include "solver.hpp"
using namespace solver;

int main()
{
    RealVariable x;
    cout << solve((x ^ 2) + 5 * x == 4) << endl;
    return 0;
}
