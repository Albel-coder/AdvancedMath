#include <iostream>

class Complex
{
private:

    double Real;
    double Imag;

public:

    Complex(double r = 0.0, double i = 0.0)
    {
        Real = r;
        Imag = i;
    }

    // Arithmetic operators
    Complex operator+(const Complex& OtherNum) const 
    {
        return Complex(Real + OtherNum.Real, Imag + OtherNum.Imag);
    }

    Complex operator-(const Complex& OtherNum) const
    {
        return Complex(Real - OtherNum.Real, Imag - OtherNum.Imag);
    }

    Complex operator*(const Complex& OtherNum) const
    {
        return Complex(Real * OtherNum.Real - Imag * OtherNum.Imag, Real * OtherNum.Imag + Imag * OtherNum.Real);
    }

    Complex operator/(const Complex& OtherNum) const
    {
        double RealResult = 0.0, ImagResult = 0.0;

        RealResult = (Real * OtherNum.Real + Imag * OtherNum.Imag) / (OtherNum.Real * OtherNum.Real + OtherNum.Imag * OtherNum.Imag);
        ImagResult = (OtherNum.Real * Imag - Real * OtherNum.Imag) / (OtherNum.Real * OtherNum.Real + OtherNum.Imag * OtherNum.Imag);

        return Complex(RealResult, ImagResult);
    }

    // Comparison Operators
    bool operator==(const Complex& OtherNum) const
    {
        if (Real == OtherNum.Real && Imag == OtherNum.Imag)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    bool operator!=(const Complex& OtherNum) const
    {
        if (Real != OtherNum.Real || Imag != OtherNum.Imag)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    // Assignment Operators
    Complex& operator=(const Complex& OtherNum)
    {
        Real = OtherNum.Real;
        Imag = OtherNum.Imag;
        return *this;
    }

    Complex& operator+=(const Complex& OtherNum)
    {
        Real += OtherNum.Real;
        Imag += OtherNum.Imag;
        return *this;
    }

    Complex& operator-=(const Complex& OtherNum)
    {
        Real -= OtherNum.Real;
        Imag -= OtherNum.Imag;
        return *this;
    }

    Complex& operator*=(const Complex& OtherNum)
    {
        double NewReal = Real * OtherNum.Real - Imag * OtherNum.Imag;
        double NewImag = Real * OtherNum.Imag + Imag * OtherNum.Real;
        Real = NewReal;
        Imag = NewImag;
        return *this;
    }

    Complex& operator/=(const Complex& OtherNum) {
        *this = *this / OtherNum;
        return *this;
    }

    // Unary operators
    Complex operator+() const 
    { 
        return *this; 
    }
    Complex operator-() const 
    { 
        return Complex(-Real, -Imag); 
    }

    // Access Methods
    double GetReal() const
    {
        return Real;
    }
    double GetImag() const
    {
        return Imag;
    }

    // Additional Methods
    double Magnitude() const 
    {
        return std::sqrt(Real * Real + Imag * Imag);
    }

    double Phase() const 
    {
        return std::atan2(Imag, Real);
    }

    Complex Conjugate() const 
    {
        return Complex(Real, -Imag);
    }

    // Friendly functions for streaming I/O
    friend std::ostream& operator<<(std::ostream& os, const Complex& c);
    friend std::istream& operator>>(std::istream& is, Complex& c);
};

// Streaming output
std::ostream& operator<<(std::ostream& Output, const Complex& C) 
{
    Output << "(" << C.Real;
    if (C.Imag >= 0)
        Output << "+" << C.Imag << "i)";
    else
        Output << C.Imag << "i)";
    return Output;
}

// Stream input
std::istream& operator>>(std::istream& Input, Complex& C) 
{
    char Plus;
    char i;
    char Paren;
    Input >> Paren >> C.Real >> Plus >> C.Imag >> i >> Paren;
    if (Plus == '-') {
        C.Imag = -C.Imag;
        Input >> i;
    }
    if (i != 'i') Input.setstate(std::ios::failbit);
    return Input;
}

int main()
{

    // Testing all operations
    Complex a(1, 2);
    Complex b(3, -4);

    std::cout << "a = " << a << "\n";
    std::cout << "b = " << b << "\n";

    std::cout << "a + b = " << a + b << "\n";
    std::cout << "a - b = " << a - b << "\n";
    std::cout << "a * b = " << a * b << "\n";
    std::cout << "a / b = " << a / b << "\n";

    Complex c = a;
    c += b;
    std::cout << "After a += b: " << c << "\n";

    c = a;
    c *= b;
    std::cout << "After a *= b: " << c << "\n";

    std::cout << "Module a: " << a.Magnitude() << "\n";
    std::cout << "Phase a: " << a.Phase() << "\n";
    std::cout << "Pairing a: " << a.Conjugate() << "\n";

    return 0;
}
