#include "PlaneTransforms.h"
#include <cmath>
#include <algorithm>
#include <sstream>

// Homothety implementation
Homothety::Homothety(const Complex& ComplexNumber, double K)
{
	Center = ComplexNumber;
	Coefficient = K;
}

Homothety::Homothety(double CenterX, double CenterY, double K)
{
	Center = Complex(CenterX, CenterY);
	Coefficient = K;
}

Complex Homothety::GetCenter() const
{
	return Center;
}

double Homothety::GetCoefficient() const
{
	return Coefficient;
}

void Homothety::SetCenter(const Complex& ComplexNumber)
{
	Center = ComplexNumber;
}

void Homothety::SetCenter(double X, double Y)
{
	Center = Complex(X, Y);
}

void Homothety::SetCoefficient(double K)
{
	Coefficient = K;
}

Complex Homothety::ApplyTo(const Complex& Point) const
{
	return Center + (Point - Center) * Coefficient;
}

Homothety Homothety::ComposeWith(const Homothety& Other) const
{
	if (Center == Other.Center)
	{
		return Homothety(Center, Coefficient * Other.Coefficient);
	}

	double NewCoefficient = Coefficient * Other.Coefficient;
	Complex NewCenter = (Center * (1 - Other.Coefficient) + Other.Center * (Other.Coefficient * (1 - Coefficient)))
		/ (1 - Coefficient * Other.Coefficient);

	return Homothety(NewCenter, NewCoefficient);
}

Homothety Homothety::getInverse() const
{
	if (Coefficient == 0.0)
	{
		return Homothety();
	}
	else
	{
		return Homothety(Center, 1.0 / Coefficient);
	}
}

Homothety Homothety::RaisedTo(int Number) const
{
	return Homothety(Center, std::pow(Coefficient, Number));
}

bool Homothety::IsExpansion() const
{
	return Coefficient > 1.0;
}

bool Homothety::IsContraction() const
{
	return Coefficient > 0.0 && Coefficient < 1.0;
}

bool Homothety::IsReflection() const
{
	return Coefficient < 0.0;
}

bool Homothety::IsIdentity() const
{
	return false;
}

double Homothety::GetScaleFactor() const
{
	return 0.0;
}

Complex Homothety::GetFixedPointers() const
{
	// Для гомотетии единственная неподвижная точка - центр (только если коэффициент не равен единице)
	if (std::abs(Coefficient - 1.0) > 1e-10)
	{
		return Center;
	}
	else // Если коэффициент равен единице, то все точки неподвижны (тождественное преобразования)
	{
		return Complex();
	}
}

std::vector<Complex> Homothety::ApplyToPolygon(const std::vector<Complex>& Polygon) const
{
	std::vector<Complex> Result;
	Result.reserve(Polygon.size());
	for (const auto& Point : Polygon)
	{
		Result.push_back(ApplyTo(Point));
	}

	return Result;
}

std::pair<Complex, double> Homothety::ApplyToCircle(const Complex& Center, double Radius) const
{
	Complex NewCenter = ApplyTo(Center);
	double NewRadius = std::abs(Coefficient) * Radius;
	return {NewCenter, NewRadius};
}

Homothety Homothety::FromScaleFactor(double Scale)
{
	return Homothety(Complex(0, 0), Scale);
}

Homothety Homothety::FromFixedPointAndImage(const Complex& FixedPoint, const Complex& ImagePoint)
{
	if (FixedPoint == ImagePoint)
	{
		return Homothety(FixedPoint, 1.0); // Identity
	}

	// Для гомотетии с центром в исправленной точкой старая функция не сработает
	// (ImagePoint = FixedPoint + K * (FixedPoint))
	// Вместо этого давайте создадим гомотетию с центром в начале координат

	double K = ImagePoint.Magnitude() / FixedPoint.Magnitude();
	return Homothety(Complex(0, 0), K);
}

std::unique_ptr<PlaneTransformation> Homothety::GetInverse() const
{
	return std::make_unique<Homothety>(getInverse());
}

bool Homothety::operator==(const Homothety& OtherNumber) const
{
	return Center == OtherNumber.Center &&
		std::abs(Coefficient - OtherNumber.Coefficient) < 1e-10;
}

bool Homothety::operator!=(const Homothety& OtherNumber) const
{
	return !(*this == OtherNumber);
}

std::string Homothety::ToString() const
{
	std::stringstream Stream;
	Stream << "Homothety(center=" << Center << ", K= " << Coefficient << ")";
	return Stream.str();
}

Homothety Homothety::FromTwoPairs(const Complex& A, const Complex& A1, const Complex& B, const Complex& B1)
{
	Complex AB = B - A;
	Complex A1B1 = B1 - A1;

	if (AB.Magnitude() < 1e-10)
	{
		return Homothety();
	}
	else
	{
		double K = A1B1.Magnitude() / AB.Magnitude();
		double CrossProduct = AB.GetReal() * A1B1.GetImag() - AB.GetImag() * A1B1.GetReal();

		if (CrossProduct < 0)
		{
			K = -K;
		}

		Complex center = (A1 - A * K) / (1 - K);
		return Homothety(center, K);
	}
}

std::ostream operator<<(std::ostream& Output, const Homothety& H)
{
	Output << H.ToString();
	return Output;
}
