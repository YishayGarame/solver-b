#include "solver.hpp"
#include <complex>
using namespace solver;

double solver::solve(RealVariable &x)
{
    double a = x.pow;
    double b = x.beforeX;
    double c = x.freeNumber;
    if (a == 0 && b == 0 && c == 0) // 0 = 0  infinty salution
    {
        throw runtime_error("There are infinity solutions for this equation");
    }
    if (a != 0 && b == 0 && c > 0) // x^2 + 5 = 0  no salution
    {
        throw runtime_error("There are no solutions for this equation");
    }
    if (a != 0 && b == 0 && c == 0) // x^2 = 0
    {
        return 0;
    }
    if (a == 0 && b != 0 && c == 0) // x = 0
    {
        return 0;
    }
    if (a == 0 & b == 0 & c != 0) // 5 = 0 no salution
    {
        throw runtime_error("There are no solutions for this equation");
    }
    //more mikrei kaze  *******************************

    double ans;
    if (a != 0)
    {

        double delta = b * b - 4 * a * c;
        if (delta < 0)
        {
            throw runtime_error("There is no salution");
        }

        else if (delta == 0)
        {
            return (-b / 2 * a);
        }
        else //delta > 0 two salutions but we will only retunr one as demanded
        {
            delta = sqrt(delta);
            ans = (-b + delta) / (2 * a);
            return ans;
        }
    }
    else //coefficient for x^2 is zero
    {
        ans = (b / (-c));
        return ans;
    }
}

double solver::solve(ComplexVariable x)
{
    return 1;
}

////////////////////////////////operators RealVarible ////////////////////////////////

//operator + //

solver::RealVariable &solver::operator+(RealVariable &x, RealVariable &y)
{
    if (x.pow == y.pow)
    {
        x.pow = x.pow + y.pow;
    }
    x.beforeX = x.beforeX + y.beforeX;
    x.freeNumber = x.freeNumber + y.freeNumber;
    return x;
}
solver::RealVariable &solver::operator+(RealVariable &x, double y)
{
    x.freeNumber = x.freeNumber + y;
    return x;
}
solver::RealVariable &solver::operator+(double y, RealVariable &x)
{
    return solver::operator+(x, y);
}

//operator - //

solver::RealVariable &solver::operator-(RealVariable &x, RealVariable &y)
{
    if (x.pow == y.pow)
    {
        x.pow = x.pow - y.pow;
    }
    x.beforeX = x.beforeX - y.beforeX;
    x.freeNumber = x.freeNumber - y.freeNumber;
    return x;
}
solver::RealVariable &solver::operator-(RealVariable &x, double y)
{
    x.freeNumber = x.freeNumber - y;
    return x;
}
solver::RealVariable &solver::operator-(double y, RealVariable &x)
{
    x.freeNumber = y - x.freeNumber;
    return x;
}

//operator * //

solver::RealVariable &solver::operator*(RealVariable &x, RealVariable &y)
{
    return x;
}
solver::RealVariable &solver::operator*(RealVariable &x, double y)
{
    return x;
}
solver::RealVariable &solver::operator*(double y, RealVariable &x)
{
    return x;
}

//operator / //

solver::RealVariable &solver::operator/(RealVariable &x, RealVariable &y)
{
    return x;
}
solver::RealVariable &solver::operator/(RealVariable &x, double y)
{
    return x;
}
solver::RealVariable &solver::operator/(double y, RealVariable &x)
{
    return x;
}

//operator ^ //

solver::RealVariable &solver::operator^(RealVariable &x, RealVariable &y)
{
    return x;
}
solver::RealVariable &solver::operator^(RealVariable &x, double y)
{
    return x;
}
solver::RealVariable &solver::operator^(double y, RealVariable &x)
{
    return x;
}

//operator == //

solver::RealVariable &solver::operator==(RealVariable &x, RealVariable &y)
{
    x.pow = x.pow - y.pow;
    x.beforeX = x.beforeX - y.beforeX;
    x.freeNumber = x.freeNumber - y.freeNumber;
    return x;
}
solver::RealVariable &solver::operator==(RealVariable &x, double y)
{
    x.freeNumber = x.freeNumber - y;
    return x;
}
solver::RealVariable &solver::operator==(double y, RealVariable &x)
{
    return x;
}

////////////////////////////////operators ComplexVariable////////////////////////////////

//operator + //

solver::ComplexVariable &solver::operator+(ComplexVariable &x, ComplexVariable &y)
{
    x.re = x.re + y.re;
    x.im = x.im = y.im;
    return x;
}
solver::ComplexVariable &solver::operator+(ComplexVariable &x, std::complex<double> y)
{
    // x.re = x.re + y.real;
    // x.im = x.im = y.imag;
    return x;
}

solver::ComplexVariable &solver::operator+(std::complex<double> y, ComplexVariable &x)
{
    // x.re = x.re + y.real;
    // x.im = x.im = y.imag;
    return x;
}
solver::ComplexVariable &solver::operator+(ComplexVariable &x, double y)
{
    x.re = x.re + y;
    return x;
}
solver::ComplexVariable &solver::operator+(double y, ComplexVariable &x)
{
    x.re = x.re + y;
    return x;
}

//operator - //

solver::ComplexVariable &solver::operator-(ComplexVariable &x, ComplexVariable &y)
{
    x.re = x.re - y.re;
    x.im = x.im = y.im;
    return x;
}
solver::ComplexVariable &operator-(ComplexVariable &x, std::complex<double> y)
{
    // x.re = x.re + y.real;
    // x.im = x.im = y.imag;
    return x;
}

solver::ComplexVariable &operator-(std::complex<double> y, ComplexVariable &x)
{
    // x.re = x.re + y.real;
    // x.im = x.im = y.imag;
    return x;
}
solver::ComplexVariable &solver::operator-(ComplexVariable &x, double y)
{
    x.re = x.re - y;
    return x;
}
solver::ComplexVariable &solver::operator-(double y, ComplexVariable &x)
{
    x.re = y - x.re;
    return x;
}

//operator * //

solver::ComplexVariable &solver::operator*(ComplexVariable &x, ComplexVariable &y)
{
    x.re = ((x.re * y.re) - (x.im * y.im));
    x.im = ((x.re * y.im) + (x.im * y.re));
    return x;
}
solver::ComplexVariable &solver::operator*(ComplexVariable &x, double y)
{
    x.re = x.re * y;
    x.im = x.im * y;
    return x;
}
solver::ComplexVariable &solver::operator*(double y, ComplexVariable &x)
{
    x.re = x.re * y;
    x.im = x.im * y;
    return x;
}

//operator / //

solver::ComplexVariable &solver::operator/(ComplexVariable &x, ComplexVariable &y)
{
    x.re = ((x.re * y.re) - (x.im * y.im));
    x.im = ((x.re * y.im) + (x.im * y.re));
    return x;
}
solver::ComplexVariable &solver::operator/(ComplexVariable &x, double y)
{
    x.re = x.re * y;
    x.im = x.im * y;
    return x;
}
solver::ComplexVariable &solver::operator/(double y, ComplexVariable &x)
{
    x.re = x.re * y;
    x.im = x.im * y;
    return x;
}

//operator ^ //

solver::ComplexVariable &solver::operator^(ComplexVariable &x, ComplexVariable &y)
{
    x.re = ((x.re * y.re) - (x.im * y.im));
    x.im = ((x.re * y.im) + (x.im * y.re));
    return x;
}
solver::ComplexVariable &solver::operator^(ComplexVariable &x, double y)
{
    x.re = x.re * y;
    x.im = x.im * y;
    return x;
}
solver::ComplexVariable &solver::operator^(double y, ComplexVariable &x)
{
    x.re = x.re * y;
    x.im = x.im * y;
    return x;
}

//operator == //

solver::ComplexVariable &solver::operator==(ComplexVariable &x, ComplexVariable &y)
{
    x.re = ((x.re * y.re) - (x.im * y.im));
    x.im = ((x.re * y.im) + (x.im * y.re));
    return x;
}
solver::ComplexVariable &solver::operator==(ComplexVariable &x, double y)
{
    x.re = x.re * y;
    x.im = x.im * y;
    return x;
}
solver::ComplexVariable &solver::operator==(double y, ComplexVariable &x)
{
    x.re = x.re * y;
    x.im = x.im * y;
    return x;
}