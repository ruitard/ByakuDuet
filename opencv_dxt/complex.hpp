// clang-format off

#pragma once

#include <cmath>

namespace cv 
{

/** @brief  A complex number class.

  The template class is similar and compatible with std::complex, however it provides slightly
  more convenient access to the real and imaginary parts using through the simple field access, as opposite
  to std::complex::real() and std::complex::imag().
*/
template<typename _Tp> class Complex
{
public:

    //! default constructor
    Complex();
    Complex( _Tp _re, _Tp _im = 0 );

    //! conjugation
    Complex conj() const;

    _Tp re, im; //< the real and the imaginary parts
};

typedef Complex<float> Complexf;
typedef Complex<double> Complexd;

template<typename _Tp> inline
Complex<_Tp>::Complex()
    : re(0), im(0) {}

template<typename _Tp> inline
Complex<_Tp>::Complex( _Tp _re, _Tp _im )
    : re(_re), im(_im) {}

template<typename _Tp> inline
Complex<_Tp> Complex<_Tp>::conj() const
{
    return Complex<_Tp>(re, -im);
}


template<typename _Tp> static inline
bool operator == (const Complex<_Tp>& a, const Complex<_Tp>& b)
{
    return a.re == b.re && a.im == b.im;
}

template<typename _Tp> static inline
bool operator != (const Complex<_Tp>& a, const Complex<_Tp>& b)
{
    return a.re != b.re || a.im != b.im;
}

template<typename _Tp> static inline
Complex<_Tp> operator + (const Complex<_Tp>& a, const Complex<_Tp>& b)
{
    return Complex<_Tp>( a.re + b.re, a.im + b.im );
}

template<typename _Tp> static inline
Complex<_Tp>& operator += (Complex<_Tp>& a, const Complex<_Tp>& b)
{
    a.re += b.re; a.im += b.im;
    return a;
}

template<typename _Tp> static inline
Complex<_Tp> operator - (const Complex<_Tp>& a, const Complex<_Tp>& b)
{
    return Complex<_Tp>( a.re - b.re, a.im - b.im );
}

template<typename _Tp> static inline
Complex<_Tp>& operator -= (Complex<_Tp>& a, const Complex<_Tp>& b)
{
    a.re -= b.re; a.im -= b.im;
    return a;
}

template<typename _Tp> static inline
Complex<_Tp> operator - (const Complex<_Tp>& a)
{
    return Complex<_Tp>(-a.re, -a.im);
}

template<typename _Tp> static inline
Complex<_Tp> operator * (const Complex<_Tp>& a, const Complex<_Tp>& b)
{
    return Complex<_Tp>( a.re*b.re - a.im*b.im, a.re*b.im + a.im*b.re );
}

template<typename _Tp> static inline
Complex<_Tp> operator * (const Complex<_Tp>& a, _Tp b)
{
    return Complex<_Tp>( a.re*b, a.im*b );
}

template<typename _Tp> static inline
Complex<_Tp> operator * (_Tp b, const Complex<_Tp>& a)
{
    return Complex<_Tp>( a.re*b, a.im*b );
}

template<typename _Tp> static inline
Complex<_Tp> operator + (const Complex<_Tp>& a, _Tp b)
{
    return Complex<_Tp>( a.re + b, a.im );
}

template<typename _Tp> static inline
Complex<_Tp> operator - (const Complex<_Tp>& a, _Tp b)
{ return Complex<_Tp>( a.re - b, a.im ); }

template<typename _Tp> static inline
Complex<_Tp> operator + (_Tp b, const Complex<_Tp>& a)
{
    return Complex<_Tp>( a.re + b, a.im );
}

template<typename _Tp> static inline
Complex<_Tp> operator - (_Tp b, const Complex<_Tp>& a)
{
    return Complex<_Tp>( b - a.re, -a.im );
}

template<typename _Tp> static inline
Complex<_Tp>& operator += (Complex<_Tp>& a, _Tp b)
{
    a.re += b; return a;
}

template<typename _Tp> static inline
Complex<_Tp>& operator -= (Complex<_Tp>& a, _Tp b)
{
    a.re -= b; return a;
}

template<typename _Tp> static inline
Complex<_Tp>& operator *= (Complex<_Tp>& a, _Tp b)
{
    a.re *= b; a.im *= b; return a;
}

template<typename _Tp> static inline
double abs(const Complex<_Tp>& a)
{
    return std::sqrt( (double)a.re*a.re + (double)a.im*a.im);
}

template<typename _Tp> static inline
Complex<_Tp> operator / (const Complex<_Tp>& a, const Complex<_Tp>& b)
{
    double t = 1./((double)b.re*b.re + (double)b.im*b.im);
    return Complex<_Tp>( (_Tp)((a.re*b.re + a.im*b.im)*t),
                        (_Tp)((-a.re*b.im + a.im*b.re)*t) );
}

template<typename _Tp> static inline
Complex<_Tp>& operator /= (Complex<_Tp>& a, const Complex<_Tp>& b)
{
    a = a / b;
    return a;
}

template<typename _Tp> static inline
Complex<_Tp> operator / (const Complex<_Tp>& a, _Tp b)
{
    _Tp t = (_Tp)1/b;
    return Complex<_Tp>( a.re*t, a.im*t );
}

template<typename _Tp> static inline
Complex<_Tp> operator / (_Tp b, const Complex<_Tp>& a)
{
    return Complex<_Tp>(b)/a;
}

template<typename _Tp> static inline
Complex<_Tp> operator /= (const Complex<_Tp>& a, _Tp b)
{
    _Tp t = (_Tp)1/b;
    a.re *= t; a.im *= t; return a;
}

}