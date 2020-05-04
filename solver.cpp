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
    if (a == 0 && b != 0 && c == 0) // 6x = 0
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
            return (-b / (2 * a));
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
        ans = ((-c) / b);
        return ans;
    }
}

////////////////////////////////operators RealVarible ////////////////////////////////

//operator + // 

solver::RealVariable &solver::operator+(RealVariable &x, RealVariable &y)
{ 
    RealVariable ans;
    ans.pow = x.pow + y.pow;
    ans.beforeX = x.beforeX + y.beforeX;
    ans.freeNumber =  x.freeNumber + y.freeNumber;
    return ans;
}
solver::RealVariable &solver::operator+(RealVariable &x, double y)
{
    RealVariable ans;
    ans.pow = x.pow;
    ans.beforeX = x.beforeX;
    ans.freeNumber =  x.freeNumber + y;
    return ans;
}
solver::RealVariable &solver::operator+(double y, RealVariable &x)
{
    RealVariable ans;
    ans.pow = x.pow;
    ans.beforeX = x.beforeX;
    ans.freeNumber =  x.freeNumber + y;
    return ans;
}
//operator - 

solver::RealVariable &solver::operator-(RealVariable &x, RealVariable &y)
{
    RealVariable ans;
    ans.pow = x.pow - y.pow;
    ans.beforeX = x.beforeX - y.beforeX;
    ans.freeNumber =  x.freeNumber - y.freeNumber;
    return ans;
}
solver::RealVariable &solver::operator-(RealVariable &x, double y)
{
    RealVariable ans;
    ans.pow = x.pow;
    ans.beforeX = x.beforeX;
    ans.freeNumber =  x.freeNumber - y;
    return ans;
}
solver::RealVariable &solver::operator-(double y, RealVariable &x)
{
    RealVariable ans;
    ans.pow = x.pow;
    ans.beforeX = x.beforeX;
    ans.freeNumber =  y - x.freeNumber;
    return ans;
}
//operator * //

solver::RealVariable &solver::operator*(RealVariable &x, RealVariable &y)
{
RealVariable ans;
   if(x.pow != 0 && y.pow != 0) 
   {
        throw runtime_error("Limit of a second defree function"); // x^4
   }
   
   if((x.pow != 0 && y.beforeX != 0 ) || (y.pow != 0 && x.beforeX != 0))
   {
        throw runtime_error("Limit of a second defree function"); x^3;
   }
    if(x.pow != 0 )
   {
       ans.pow = x.pow * y.freeNumber;
       ans.beforeX = x.beforeX * y.freeNumber;
       ans.freeNumber = x.freeNumber * y.freeNumber;
   }
    else if(y.pow != 0)
   {
       ans.pow = y.pow * x.freeNumber;
       ans.beforeX = y.beforeX * x.freeNumber;
       ans.freeNumber = y.freeNumber * x.freeNumber;
   }
   else 
   {
       ans.pow = x.beforeX * y.beforeX;
       ans.beforeX = x.beforeX * y.freeNumber + x.freeNumber * y.beforeX;
       ans.freeNumber = x.freeNumber * y.freeNumber;
   }
   return ans;
}
solver::RealVariable &solver::operator*(RealVariable &x, double y)
{
    RealVariable ans;
    ans.pow = x.pow * y;
    ans.beforeX = x.beforeX * y;
    ans.freeNumber = x.freeNumber * y;
    return ans;
}
solver::RealVariable &solver::operator*(double y, RealVariable &x)
{
    RealVariable ans;
    ans.pow = x.pow * y;
    ans.beforeX = x.beforeX * y;
    ans.freeNumber = x.freeNumber * y;
    return ans;

}
// solver::RealVariable &solver::operator/(RealVariable &x, RealVariable &y)
// {
//     RealVariable ans;

//     return ans;
// }
solver::RealVariable &solver::operator/(RealVariable &x, double y)
{
    RealVariable ans;
    if (y == 0){
        throw runtime_error("cant divid by zero"); 
    }
    ans.pow = x.pow / y;
    ans.beforeX = x.beforeX / y;
    ans.freeNumber = x.freeNumber / y;  

    return ans;
}
solver::RealVariable &solver::operator/(double y, RealVariable &x)
{
    RealVariable ans;
    if (x.pow == 0 || x.beforeX == 0 || x.freeNumber == 0){
        throw runtime_error("cant divid by zero"); 
    }
    ans.pow = y / x.pow;
    ans.beforeX = y / x.beforeX;
    ans.freeNumber = y / x.freeNumber;  
    return ans;
}

//operator ^ //

// solver::RealVariable &solver::operator^(RealVariable &x, RealVariable &y)
// {
//     return x;
// }
solver::RealVariable &solver::operator^(RealVariable &x, double y)
{
        RealVariable ans;
        if (y > 2 || y < 0 || ((y == 2 && x.pow != 0)))
        {
            throw runtime_error("limit degree of 0-2");
        }

        else if (y == 0)      // y = 0
        {
            ans.pow == 0;
            ans.beforeX = 0;
            ans.freeNumber = 1;
        }
        else if (y == 1)            // y = 1
        {        
            ans.pow = x.pow;
            ans.beforeX = x.beforeX;
            ans.freeNumber = x.freeNumber;
        } 
        else if (y == 2 && x.beforeX == 0)        // 5^2
        {
            ans.pow = 0;
            ans.beforeX = 0;
            ans.freeNumber = x.freeNumber * x.freeNumber;
        }
          else if (y == 2 && x.freeNumber == 0)      // x^2
        {
            ans.pow = x.beforeX * x.beforeX;
            ans.beforeX = 0;
            ans.freeNumber = 0;
        } 
        else                                    // (x + 3) ^ 2
        {
            ans.pow = x.beforeX * x.beforeX;
            ans.beforeX = x.beforeX * x.freeNumber;
            ans.freeNumber = x.freeNumber * x.freeNumber;
        }

        return ans;

}
// solver::RealVariable &solver::operator^(double y, RealVariable &x)
// {
//     return x;
// }

//operator == //

solver::RealVariable &solver::operator==(RealVariable &x, RealVariable &y)
{
    RealVariable ans;
    ans.pow = x.pow - y.pow;
    ans.beforeX = x.beforeX - y.beforeX;
    ans.freeNumber = x.freeNumber - y.freeNumber;
    return ans;
}
solver::RealVariable &solver::operator==(RealVariable &x, double y)
{
    RealVariable ans;
    ans.freeNumber = x.freeNumber - y;
    return ans;
}
solver::RealVariable &solver::operator==(double y, RealVariable &x)
{
    RealVariable ans;
    ans.pow = - x.pow;
    ans.beforeX = - x.beforeX;
    ans.freeNumber =  y - x.freeNumber;
    return ans;
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

double solver::solve(ComplexVariable x)
{
    return 1;
}