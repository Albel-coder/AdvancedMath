#pragma once
#include "Complex.h"
#include <stdexcept>
#include <string>
#include <vector>
#include <memory>

// Ѕазовый класс дл€ всех преобразований плоскости
class PlaneTransformation
{
public:
	virtual ~PlaneTransformation() = default;
	virtual Complex ApplyTo(const Complex& Point) const = 0;
	virtual std::string ToString() const = 0;
	virtual bool IsIdentity() const = 0;
	virtual std::unique_ptr<PlaneTransformation> GetInverse() const = 0;
	virtual std::unique_ptr<PlaneTransformation> ComposeWith(const PlaneTransformation& Other) const = 0;
};

class Homothety : public PlaneTransformation
{
private:
	Complex Center;
	double Coefficient;

public:
	Homothety(const Complex& ComplexNumber = Complex(0.0, 0.0), double K = 1.0);
	Homothety(double CenterX, double CenterY, double K);

	// Basic getters/setters
	Complex GetCenter() const;
	double GetCoefficient() const;
	void SetCenter(const Complex& ComplexNumber);
	void SetCenter(double X, double Y);
	void SetCoefficient(double K);

	// Transformation methods
	Complex ApplyTo(const Complex& Point) const override;
	Homothety ComposeWith(const Homothety& Other) const;
	Homothety getInverse() const;
	Homothety RaisedTo(int Number) const;

	// New advanced methods
	bool IsExpansion() const;      // k > 1
	bool IsContraction() const;    // 0 < k < 1
	bool IsReflection() const;     // k < 0
	bool IsIdentity() const override;

	double GetScaleFactor() const;
	Complex GetFixedPointers() const; // ¬озвращает неподвижные точки

	// Apply to geometric objects
	std::vector<Complex> ApplyToPolygon(const std::vector<Complex>& Polygon) const;
	std::pair<Complex, double> ApplyToCircle(const Complex& Center, double Radius) const;

	// Static factory methods
	static Homothety FromScaleFactor(double Scale);
	static Homothety FromFixedPointAndImage(const Complex& FixedPoint, const Complex& ImagePoint);

	// PlaneTransformation interface
	std::unique_ptr<PlaneTransformation> GetInverse() const override;
	std::unique_ptr<PlaneTransformation> ComposeWith(const PlaneTransformation& Other) const override;

	// Operators
	bool operator==(const Homothety& OtherNumber) const;
	bool operator!=(const Homothety& OtherNumber) const;
	Homothety operator*(const Homothety& OtherNumber) const;

	std::string ToString() const override;

	static Homothety FromTwoPairs(const Complex& A, const Complex& A1,
		const Complex& B, const Complex& B1);
};

std::ostream operator<<(std::ostream& Output, const Homothety& H);