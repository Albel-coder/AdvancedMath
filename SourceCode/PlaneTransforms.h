#pragma once
#include "Complex.h"
#include <stdexcept>
#include <string>
#include <vector>
#include <memory>
#include <tuple>

// Base class for all plane transformations
class PlaneTransformation
{
public:
	virtual ~PlaneTransformation() = default;
	virtual Complex applyTo(const Complex& point) const = 0;
	virtual std::string toString() const = 0;
	virtual bool isIdentity() const = 0;
	virtual std::unique_ptr<PlaneTransformation> getInverse() const = 0;
	virtual std::unique_ptr<PlaneTransformation> composeWith(const PlaneTransformation& other) const = 0;
};

class Homothety : public PlaneTransformation
{
private:
	Complex center;
	double coefficient;

public:
	Homothety(const Complex& complexNumber = Complex(0.0, 0.0), double k = 1.0);
	Homothety(double centerX, double centerY, double k);

	// Basic getters/setters
	Complex getCenter() const;
	double getCoefficient() const;
	void setCenter(const Complex& complexNumber);
	void setCenter(double x, double y);
	void setCoefficient(double k);

	// Transformation methods
	Complex applyTo(const Complex& point) const override;
	Homothety composeWith(const Homothety& other) const;
	Homothety getInverse();
	Homothety raisedTo(int number) const;

	// New advanced methods
	bool isExpansion() const;      // k > 1
	bool isContraction() const;    // 0 < k < 1
	bool isReflection() const;     // k < 0
	bool isIdentity() const override;

	double getScaleFactor() const;
	Complex getFixedPointers() const; // Returns fixed points

	// Apply to geometric objects
	std::vector<Complex> applyToPolygon(const std::vector<Complex>& polygon) const;
	std::pair<Complex, double> applyToCircle(const Complex& center, double radius) const;

	// Static factory methods
	static Homothety fromScaleFactor(double scale);
	static Homothety fromFixedPointAndImage(const Complex& fixedPoint, const Complex& imagePoint);

	// PlaneTransformation interface
	std::unique_ptr<PlaneTransformation> getInverse() const override;
	std::unique_ptr<PlaneTransformation> composeWith(const PlaneTransformation& other) const override;

	// Operators
	bool operator==(const Homothety& otherNumber) const;
	bool operator!=(const Homothety& otherNumber) const;
	Homothety operator*(const Homothety& otherNumber) const;

	std::string toString() const override;

	static Homothety fromTwoPairs(const Complex& a, const Complex& a1,
		const Complex& b, const Complex& b1);
};

std::ostream& operator<<(std::ostream& output, const Homothety& h);

class AffineTransform : public PlaneTransformation
{
private:
	double Matrix11;
	double Matrix12;
	double Matrix13;
	double Matrix21;
	double Matrix22;
	double Matrix23;

public:
	AffineTransform(double matrix11 = 1.0, double matrix12 = 0.0, double matrix13 = 0.0,
		double matrix21 = 0.0, double matrix22 = 1.0, double matrix23 = 0.0);

	// Static factory methods
	static AffineTransform Identity();
	static AffineTransform Translation(double tx, double ty);
	static AffineTransform Translation(const Complex& Vector);
	static AffineTransform Scaling(double sx, double sy);
	static AffineTransform Scaling(double Scale);
	static AffineTransform Scaling(const Complex& Center, double sx, double sy);
	static AffineTransform Rotation(double AngleRadius);
	static AffineTransform Rotation(double AngleRadius, const Complex& Center);
	static AffineTransform RotationDegrees(double AngleDegrees);
	static AffineTransform RotationDegrees(double AngleDegrees, const Complex& Center);
	static AffineTransform Shear(double shx, double shy);
	static AffineTransform ReflectionOverLine(const Complex& FirstPoint, const Complex& SecondPoint);
	static AffineTransform ReflectionOverXAxis();
	static AffineTransform ReflectionOverYAxis();
	static AffineTransform ProjectionOntoLine(const Complex& FirstPoint, const Complex& SecondPoint);

	// Advanced static constructors
	static AffineTransform FromTriangleMapping(const Complex& srcA, const Complex& srcB, const Complex& srcC,
		const Complex& dstA, const Complex& dstB, const Complex& dstC);
	
	static AffineTransform FromBasisVectors(const Complex& FirstE, const Complex& SecondE, const Complex& Origin = Complex(0.0, 0.0));

	// Transformation methods
	Complex applyTo(const Complex& Point) const override;
	AffineTransform ComposeWith(const AffineTransform& Other) const;
	AffineTransform operator*(const AffineTransform& Other) const;
	AffineTransform Inverse() const;

	// New advanced methods
	bool isIdentity() const override;
	bool IsIsometry() const;
	bool IsSimilarity() const;    // Keeping the angles
	bool IsEquiareal() const;     // We preserve space
	bool IsDirect() const;        // Maintaining orientation
	bool IsInvolutory() const;    // T^2 = Identity

	// Geometric properties
	Complex GetTranslation() const;
	double GetRotationAngle() const;
	Complex GetScaleFactors() const;
	Complex GetShearFactors() const;
	double GetAreaScale() const;     // Area change factor

	// Fixed points and invariant lines
	std::vector<Complex> GetFixedPoints() const;
	bool HasFixedPoint(const Complex& Point) const;

	// Eigen analysis
	std::vector<Complex> GetEigenvalues() const;
	std::vector<Complex> GetEigenvectors() const;

	// Apply to geometric objects
	std::vector<Complex> ApplyToPolygon(const std::vector<Complex>& Polygon) const;
	std::pair<Complex, double> ApplyToCircle(const Complex& Center, double Radius) const;
	Complex ApplyToVector(const Complex& Vector) const; // No broadcast

	// Decomposition methods
	void Decompose(Complex& Translation, Complex& Scale, double& RotationAngle) const;
	std::tuple<AffineTransform, AffineTransform, AffineTransform> DecomposeLDU() const; // LDU decomposition
	std::tuple<AffineTransform, AffineTransform, AffineTransform> DecomposeQR() const; // QR decomposition

	// Matrix operations
	AffineTransform Transpose() const;
	double Trace() const;
	double Determinant() const;
	AffineTransform Adjugate() const;

	// Power and exponential
	AffineTransform Power(int Number) const;
	AffineTransform Exponential() const;

	// Interpolation
	static AffineTransform Lerp(const AffineTransform& A, const AffineTransform& B, double t);

	// PlaneTransformation interface
	std::unique_ptr<PlaneTransformation> getInverse() const override;
	std::unique_ptr<PlaneTransformation> composeWith(const PlaneTransformation& Other) const override;

	// Getters
	void GetMatrix(double& Matrix11, double& Matrix12, double& Matrix13,
		double& Matrix21, double& Matrix22, double& Matrix23) const;
	void GetLinearPart(double& Matrix11, double& Matrix12, double& Matrix21, double& Matrix22) const;
	void GetTranslation(double& tx, double& ty) const;

	std::string toString() const override;
	std::string ToMatrixString() const;

	bool operator==(const AffineTransform& Other) const;
	bool operator!=(const AffineTransform& Other) const;

	};

	class CompositeTransformation : public PlaneTransformation {
	private:
		std::vector<std::unique_ptr<PlaneTransformation>> transforms;
	public:
		void AddTransformation(std::unique_ptr<PlaneTransformation> t);
		Complex applyTo(const Complex& Point) const override;
		std::string toString() const override;
		bool isIdentity() const override;
		std::unique_ptr<PlaneTransformation> getInverse() const override;
		std::unique_ptr<PlaneTransformation> composeWith(const PlaneTransformation& Other) const override;
	};

	std::ostream& operator<<(std::ostream& Output, const AffineTransform& Transform);