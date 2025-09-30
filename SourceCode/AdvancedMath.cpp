#include <iostream>
#include <cmath>
#include <numbers>
#include <ostream>
#include <string>
#include <vector>

const double PI = 3.1415926535897932;

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

    // Convert to trigonometric form
    ComplexTrigonometric toTrigonometric() const 
    {
        return ComplexTrigonometric::FromAlgebraic(Real, Imag);
    }

    // Creation from trigonometric form
    static Complex FromTrigonometric(const ComplexTrigonometric& CT) 
    {
        return Complex(CT.toReal(), CT.toImaginary());
    }

    // Operations with trigonometric numbers
    Complex operator*(const ComplexTrigonometric& OtherNumber) const 
    {
        return *this * FromTrigonometric(OtherNumber);
    }

    Complex operator/(const ComplexTrigonometric& OtherNumber) const 
    {
        return *this / FromTrigonometric(OtherNumber);
    }

    // Friendly functions for streaming I/O
    friend std::ostream& operator<<(std::ostream& Output, const Complex& C);
    friend std::istream& operator>>(std::istream& Input, Complex& C);
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
    if (Plus == '-') 
    {
        C.Imag = -C.Imag;
        Input >> i;
    }
    if (i != 'i') Input.setstate(std::ios::failbit);
    return Input;
}

class ComplexTrigonometric {
private:

    double Magnitude;  // modulus (r)
    double Angle;      // argument (θ) in radians


public:

    // Constructors
    ComplexTrigonometric(double mag = 0.0, double ang = 0.0)
        : Magnitude(mag), Angle(ang) { }

    // Creation from polar coordinates
    static ComplexTrigonometric FromPolar(double Mag, double Ang) 
    {
        return ComplexTrigonometric(Mag, Ang);
    }

    // Creating from polar coordinates with an angle in degrees
    static ComplexTrigonometric fromPolarDegrees(double Mag, double Degrees) 
    {
        return ComplexTrigonometric(Mag, Degrees * PI / 180.0);
    }

    // Creation from algebraic form
    static ComplexTrigonometric FromAlgebraic(double Real, double Imag);

    // Creation from an ordinary complex number
    static ComplexTrigonometric FromComplex(const Complex& C);

    // Getters
    double getMagnitude() const 
    { 
        return Magnitude; 
    }
    double getAngle() const 
    { 
        return Angle; 
    }
    double getAngleDegrees() const 
    { 
        return Angle * 180.0 / PI; 
    }

    // Setters
    void setMagnitude(double Mag) 
    { 
        Magnitude = Mag; 
    }
    void setAngle(double Ang) 
    { 
        Angle = Ang; 
    }
    void setAngleDegrees(double Degrees) {
        Angle = Degrees * PI / 180.0;
    }

    // Arithmetic operations in trigonometric form
    ComplexTrigonometric operator*(const ComplexTrigonometric& Other) const 
    {
        // Multiplication: modules are multiplied, angles are added
        return ComplexTrigonometric(
            Magnitude * Other.Magnitude,
            Angle + Other.Angle
        );
    }

    ComplexTrigonometric operator/(const ComplexTrigonometric& Other) const 
    {
        // Division: modules are divided, angles are subtracted
        if (Other.Magnitude == 0.0) {
            throw std::runtime_error("Division by zero");
        }
        return ComplexTrigonometric(
            Magnitude / Other.Magnitude,
            Angle - Other.Angle
        );
    }

    ComplexTrigonometric operator^(int32_t Power) const 
    {
        // Exponentiation: modulus to the power, angle multiplied by the power
        return ComplexTrigonometric(
            std::pow(Magnitude, Power),
            Angle * Power
        );
    }

    ComplexTrigonometric operator^(double Power) const 
    {
        // Raising to a real power
        return ComplexTrigonometric(
            std::pow(Magnitude, Power),
            Angle * Power
        );
    }

    // Root extraction
    std::vector<ComplexTrigonometric> roots(int Number) const 
    {
        if (Number <= 0) 
        {
            throw std::runtime_error("Root degree must be positive");
        }

        std::vector<ComplexTrigonometric> result;
        double root_magnitude = std::pow(Magnitude, 1.0 / Number);

        for (int k = 0; k < Number; ++k) 
        {
            double root_angle = (Angle + 2 * 3.1415 * k) / Number;
            result.emplace_back(root_magnitude, root_angle);
        }

        return result;
    }

    // Conversion to algebraic form
    double toReal() const 
    {
        return Magnitude * std::cos(Angle);
    }

    double toImaginary() const 
    {
        return Magnitude * std::sin(Angle);
    }

    // Obtaining the conjugate number
    ComplexTrigonometric conjugate() const 
    {
        return ComplexTrigonometric(Magnitude, -Angle);
    }

    // Normalize an angle to the range [-π, π]
    ComplexTrigonometric normalized() const 
    {
        double normalized_angle = std::fmod(Angle, 2 * PI);
        if (normalized_angle > PI) 
        {
            normalized_angle -= 2 * PI;
        }
        else if (normalized_angle <= -PI) 
        {
            normalized_angle += 2 * PI;
        }
        return ComplexTrigonometric(Magnitude, normalized_angle);
    }

    // String representation
    std::string toString() const 
    {
        return std::to_string(Magnitude) + "(cos(" +
            std::to_string(Angle) + ") + i*sin(" +
            std::to_string(Angle) + "))";
    }

    std::string toStringDegrees() const 
    {
        return std::to_string(Magnitude) + "(cos(" +
            std::to_string(getAngleDegrees()) + "°) + i*sin(" +
            std::to_string(getAngleDegrees()) + "°))";
    }

    // Comparison Operators
    bool operator==(const ComplexTrigonometric& other) const 
    {
        auto norm1 = normalized();
        auto norm2 = other.normalized();
        return std::abs(norm1.Magnitude - norm2.Magnitude) < 1e-10 &&
            std::abs(norm1.Angle - norm2.Angle) < 1e-10;
    }

    bool operator!=(const ComplexTrigonometric& other) const 
    {
        return !(*this == other);
    }
};

// Inference operator
std::ostream& operator<<(std::ostream& Output, const ComplexTrigonometric& CT) 
{
    Output << CT.toString();
    return Output;
}

int main()
{
    // Creating trigonometric complex numbers
    ComplexTrigonometric ct1 = ComplexTrigonometric::FromPolar(2.0, PI / 3);
    ComplexTrigonometric ct2 = ComplexTrigonometric::fromPolarDegrees(3.0, 45.0);

    std::cout << "ct1 = " << ct1 << "\n";
    std::cout << "ct2 = " << ct2.toStringDegrees() << "\n";

    // Operations in trigonometric form
    ComplexTrigonometric product = ct1 * ct2;
    ComplexTrigonometric power = ct1 ^ 3;

    std::cout << "ct1 * ct2 = " << product << "\n";
    std::cout << "ct1^3 = " << power << "\n";

    // Root extraction
    auto cubeRoots = ct1.roots(3);
    std::cout << "Cube roots of ct1:" << "\n";
    for (const auto& root : cubeRoots) {
        std::cout << "  " << root << "\n";
    }

    // Conversion between forms
    Complex algebraic = Complex::FromTrigonometric(ct1);
    ComplexTrigonometric backToTrig = algebraic.toTrigonometric();

    std::cout << "Algebraic form: " << algebraic << "\n";
    std::cout << "Back to trigonometric: " << backToTrig << "\n";

    return 0;
}
