#include "Complex.h"

const double PI = 3.1415926535897932;

Complex::Complex(double r = 0.0, double i = 0.0)
{
    Real = r;
    Imag = i;
}

// Arithmetic operators
Complex Complex::operator+(const Complex& OtherNumber) const
{
    return Complex(Real + OtherNumber.Real, Imag + OtherNumber.Imag);
}

Complex Complex::operator-(const Complex& OtherNumber) const
{
    return Complex(Real - OtherNumber.Real, Imag - OtherNumber.Imag);
}

Complex Complex::operator*(const Complex& OtherNumber) const
{
    return Complex(Real * OtherNumber.Real - Imag * OtherNumber.Imag,
        Real * OtherNumber.Imag + Imag * OtherNumber.Real);
}

Complex Complex::operator/(const Complex& OtherNumber) const
{
    double RealResult = 0.0;
    double ImagResult = 0.0;

    if ((OtherNumber.Real * OtherNumber.Real + OtherNumber.Imag * OtherNumber.Imag) != 0)
    {
        RealResult = (Real * OtherNumber.Real + Imag * OtherNumber.Imag) / (OtherNumber.Real * OtherNumber.Real + OtherNumber.Imag * OtherNumber.Imag);
        ImagResult = (OtherNumber.Real * Imag - Real * OtherNumber.Imag) / (OtherNumber.Real * OtherNumber.Real + OtherNumber.Imag * OtherNumber.Imag);
    }

    return Complex(RealResult, ImagResult);
}

// Comparison Operators
bool Complex::operator==(const Complex& OtherNumber) const
{
    if (Real == OtherNumber.Real && Imag == OtherNumber.Imag)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool Complex::operator!=(const Complex& OtherNumber) const
{
    if (Real != OtherNumber.Real || Imag != OtherNumber.Imag)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool Complex::operator>(const Complex& OtherNumber) const
{
    return this->MagnitudeSqr() > OtherNumber.MagnitudeSqr();
}

bool Complex::operator<(const Complex& OtherNumber) const
{
    return this->MagnitudeSqr() < OtherNumber.MagnitudeSqr();
}

// Assignment Operators
Complex& Complex::operator=(const Complex& OtherNumber)
{
    if (this != &OtherNumber)
    {
        Real = OtherNumber.Real;
        Imag = OtherNumber.Imag;
    }
    return *this;
}

Complex& Complex::operator+=(const Complex& OtherNumber)
{
    Real += OtherNumber.Real;
    Imag += OtherNumber.Imag;
    return *this;
}

Complex& Complex::operator-=(const Complex& OtherNumber)
{
    Real -= OtherNumber.Real;
    Imag -= OtherNumber.Imag;
    return *this;
}

Complex& Complex::operator*=(const Complex& OtherNumber)
{
    double NewReal = Real * OtherNumber.Real - Imag * OtherNumber.Imag;
    double NewImag = Real * OtherNumber.Imag + Imag * OtherNumber.Real;
    Real = NewReal;
    Imag = NewImag;
    return *this;
}

Complex& Complex::operator/=(const Complex& OtherNumber)
{
    if (OtherNumber != 0)
    {
        *this = *this / OtherNumber;
    }

    return *this;
}

// Unary operators
Complex Complex::operator+()
{
    return *this;
}

Complex Complex::operator-()
{
    return Complex(-Real, -Imag);
}

// Access Methods
double Complex::GetReal() const
{
    return Real;
}

double Complex::GetImag() const
{
    return Imag;
}

void Complex::SetReal(double r)
{
    Real = r;
}

void Complex::SetImage(double i)
{
    Imag = i;
}

// Basic mathematical operations
double Complex::Magnitude() const
{
    return std::sqrt(Real * Real + Imag * Imag);
}

double Complex::Phase() const
{
    return std::atan2(Imag, Real);
}

Complex Complex::Conjugate() const
{
    return Complex(Real, -Imag);
}

double Complex::MagnitudeSqr() const
{
    return Real * Real + Imag * Imag;
}

// External operators for working with doubles
Complex operator+(double ihs, const Complex rhs)
{
    return Complex(ihs) + rhs;
}

Complex operator-(double ihs, const Complex rhs)
{
    return Complex(ihs) - rhs;
}

Complex operator*(double ihs, const Complex rhs)
{
    return Complex(ihs) * rhs;
}

Complex operator/(double ihs, const Complex rhs)
{
    if (rhs != 0)
    {
        return Complex(ihs) / rhs;
    }
    else
    {
        return Complex(0.0, 0.0);
    }
}
