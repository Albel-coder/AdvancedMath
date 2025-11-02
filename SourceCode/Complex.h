#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <string>

class Complex
{
private:
    double Real;
    double Imag;

public:
    Complex(double r = 0.0, double i = 0.0);

    // Arithmetic operators
    Complex operator+(const Complex& OtherNumber) const;
    Complex operator-(const Complex& OtherNumber) const;
    Complex operator*(const Complex& OtherNumber) const;
    Complex operator/(const Complex& OtherNumber) const;

    // Comparison Operators
    bool operator==(const Complex& OtherNumber) const;
    bool operator!=(const Complex& OtherNumber) const;
    bool operator>(const Complex& OtherNumber) const;
    bool operator<(const Complex& OtherNumber) const;

    // Assignment Operators
    Complex& operator=(const Complex& OtherNumber);
    Complex& operator+=(const Complex& OtherNumber);
    Complex& operator-=(const Complex& OtherNumber);
    Complex& operator*=(const Complex& OtherNumber);
    Complex& operator/=(const Complex& OtherNumber);

    Complex operator+();
    Complex operator-();

    // Access Methods
    double GetReal() const;
    double GetImag() const;
    void SetReal(double r);
    void SetImage(double i);

    // Additional Methods
    double Magnitude() const;
    double Phase() const;
    Complex Conjugate() const;
    double MagnitudeSqr() const;  
};

// Streaming output
std::ostream& operator<<(std::ostream& Output, const Complex& ComplexNumber);