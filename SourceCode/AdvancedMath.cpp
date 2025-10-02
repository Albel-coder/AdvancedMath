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
        : Magnitude(mag), Angle(ang) {
    }

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

class Homothety 
{
private:

    Complex Center;
    double Coefficient;

public:
    // Constructors
    Homothety(const Complex& C = Complex(0.0, 0.0), double K = 1.0)
        : Center(C), Coefficient(K) {}

    Homothety(double CenterX, double CenterY, double k)
        : Center(CenterX, CenterY), Coefficient(k) {}

    // Getters
    Complex GetCenter() const 
    {
        return Center; 
    }
    double GetCoefficient() const 
    {
        return Coefficient; 
    }

    // Setters
    void SetCenter(const Complex& c) 
    { 
        Center = c; 
    }
    void SetCenter(double x, double y) 
    { 
        Center = Complex(x, y); 
    }
    void SetCoefficient(double k) 
    { 
        Coefficient = k; 
    }

    // Apply homothety to a point
    Complex ApplyTo(const Complex& point) const 
    {
        return Center + (point - Center) * Coefficient;
    }

    // Apply homothety to a trigonometric complex number
    ComplexTrigonometric ApplyTo(const ComplexTrigonometric& point) const 
    {
        Complex algebraicPoint = Complex::FromTrigonometric(point);
        Complex result = ApplyTo(algebraicPoint);
        return result.toTrigonometric();
    }

    // Composition of homotheties (now named more clearly)
    Homothety ComposeWith(const Homothety& other) const 
    {
        // H2 ∘ H1 = Homothety with new center and coefficient
        double newCoefficient = Coefficient * other.Coefficient;
        Complex newCenter = Center + (other.Center - Center) * (Coefficient / (1 - Coefficient * other.Coefficient));

        return Homothety(newCenter, newCoefficient);
    }

    // Inverse homothety
    Homothety getInverse() const 
    {
        if (Coefficient == 0.0) 
        {
            throw std::runtime_error("Homothety with coefficient 0 is not invertible");
        }
        return Homothety(Center, 1.0 / Coefficient);
    }

    // Power of homothety
    Homothety raisedTo(int Number) const 
    {
        return Homothety(Center, std::pow(Coefficient, Number));
    }

    // Comparison operators
    bool operator==(const Homothety& OtherNumber) const {
        return Center == OtherNumber.Center &&
            std::abs(Coefficient - OtherNumber.Coefficient) < 1e-10;
    }

    bool operator!=(const Homothety& OtherNumber) const 
    {
        return !(*this == OtherNumber);
    }

    // String representation
    std::string toString() const 
    {
        std::string realStr = std::to_string(Center.GetReal());
        std::string imagStr = std::to_string(Center.GetImag());
        std::string kStr = std::to_string(Coefficient);

        // Remove trailing zeros
        realStr = realStr.substr(0, realStr.find('.') + 3);
        imagStr = imagStr.substr(0, imagStr.find('.') + 3);
        kStr = kStr.substr(0, kStr.find('.') + 3);

        return "Homothety(center=" + realStr +
            (Center.GetImag() >= 0 ? "+" : "") + imagStr + "i, k=" + kStr + ")";
    }

    // Create homothety from two pairs of points (A->A', B->B')
    static Homothety fromTwoPairs(const Complex& A, const Complex& A1,
        const Complex& B, const Complex& B1) 
    {
        Complex AB = B - A;
        Complex A1B1 = B1 - A1;

        if (AB.Magnitude() < 1e-10) 
        {
            throw std::runtime_error("Points A and B are too close");
        }

        // k is the ratio of distances
        double k = A1B1.Magnitude() / AB.Magnitude();

        // Determine sign based on orientation
        double crossProduct = AB.GetReal() * A1B1.GetImag() - AB.GetImag() * A1B1.GetReal();
        if (crossProduct < 0) k = -k;

        // Calculate center from the equation: A1 = center + k*(A - center)
        // => center = (A1 - k*A) / (1 - k)
        Complex center = (A1 - A * k) / (1 - k);

        return Homothety(center, k);
    }
};

// Output operator
std::ostream& operator<<(std::ostream& os, const Homothety& h) {
    os << h.toString();
    return os;
}

class AffineTransform 
{
private:

    double Matrix11, Matrix12, Matrix13;  // 2x3 transformation matrix
    double Matrix21, Matrix22, Matrix23;  // [Matrix11 Matrix12 Matrix13]
    // [Matrix21 Matrix22 Matrix23]

public:

    // Constructors
    AffineTransform(double a11 = 1.0, double a12 = 0.0, double a13 = 0.0,
        double a21 = 0.0, double a22 = 1.0, double a23 = 0.0)
        : Matrix11(a11), Matrix12(a12), Matrix13(a13), Matrix21(a21), Matrix22(a22), Matrix23(a23) {
    }

    // Unit transformation
    static AffineTransform Identity() 
    {
        return AffineTransform(1.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    }

    // Creating a Shear Transformation
    static AffineTransform Translation(double tx, double ty) 
    {
        return AffineTransform(1.0, 0.0, tx, 0.0, 1.0, ty);
    }

    static AffineTransform Translation(const Complex& vector) 
    {
        return Translation(vector.GetReal(), vector.GetImag());
    }

    // Creating a Scaling Transform
    static AffineTransform Scaling(double sx, double sy) 
    {
        return AffineTransform(sx, 0.0, 0.0, 0.0, sy, 0.0);
    }

    static AffineTransform Scaling(double scale) 
    {
        return Scaling(scale, scale);
    }

    static AffineTransform Scaling(const Complex& center, double sx, double sy) 
    {
        return Translation(center) * Scaling(sx, sy) * Translation(-center);
    }

    // Creating a Rotation Transform
    static AffineTransform Rotation(double angle_radians) 
    {
        double cos_a = std::cos(angle_radians);
        double sin_a = std::sin(angle_radians);
        return AffineTransform(cos_a, -sin_a, 0.0, sin_a, cos_a, 0.0);
    }

    static AffineTransform Rotation(double angle_radians, const Complex& center) 
    {
        return Translation(center) * Rotation(angle_radians) * Translation(-center);
    }

    static AffineTransform RotationDegrees(double angle_degrees) 
    {
        return Rotation(angle_degrees * PI / 180.0);
    }

    static AffineTransform RotationDegrees(double angle_degrees, const Complex& center) 
    {
        return Rotation(angle_degrees * PI / 180.0, center);
    }

    // Creating a shear transformation
    static AffineTransform Shear(double shx, double shy) 
    {
        return AffineTransform(1.0, shx, 0.0, shy, 1.0, 0.0);
    }

    // Applying a transformation to a point
    Complex ApplyTo(const Complex& point) const 
    {
        double x = point.GetReal();
        double y = point.GetImag();
        double new_x = Matrix11 * x + Matrix12 * y + Matrix13;
        double new_y = Matrix21 * x + Matrix22 * y + Matrix23;
        return Complex(new_x, new_y);
    }

    // Applying a transformation to a trigonometric number
    ComplexTrigonometric ApplyTo(const ComplexTrigonometric& point) const 
    {
        Complex algebraic_point = Complex::FromTrigonometric(point);
        Complex result = ApplyTo(algebraic_point);
        return result.toTrigonometric();
    }

    // Composition of transformations (matrix multiplication)
    AffineTransform ComposeWith(const AffineTransform& other) const 
    {
        return AffineTransform(
            Matrix11 * other.Matrix11 + Matrix12 * other.Matrix21,      // New Matrix11
            Matrix11 * other.Matrix12 + Matrix12 * other.Matrix22,      // New Matrix12
            Matrix11 * other.Matrix13 + Matrix12 * other.Matrix23 + Matrix13, // New Matrix13
            Matrix21 * other.Matrix11 + Matrix22 * other.Matrix21,      // New Matrix21
            Matrix21 * other.Matrix12 + Matrix22 * other.Matrix22,      // New Matrix22
            Matrix21 * other.Matrix13 + Matrix22 * other.Matrix23 + Matrix23 // New Matrix23
        );
    }

    // Composition operator
    AffineTransform operator*(const AffineTransform& other) const 
    {
        return ComposeWith(other);
    }

    // Reverse conversion
    AffineTransform Inverse() const
    {
        double det = Matrix11 * Matrix22 - Matrix12 * Matrix21;

        if (std::abs(det) < 1e-10) {
            throw std::runtime_error("Affine transform is not invertible");
        }

        double inv_det = 1.0 / det;

            return AffineTransform(
                Matrix22 * inv_det,                          // New Matrix11
                -Matrix12 * inv_det,                         // New Matrix12
                (Matrix12 * Matrix23 - Matrix22 * Matrix13) * inv_det,      // New Matrix13
                -Matrix21 * inv_det,                         // New Matrix21
                Matrix11 * inv_det,                          // New Matrix22
                (Matrix21 * Matrix13 - Matrix11 * Matrix23) * inv_det       // New Matrix23
            );
    }

    // Determinant of the transformation matrix
    double Determinant() const 
    {
        return Matrix11 * Matrix22 - Matrix12 * Matrix21;
    }

    // Checking for orientation preservation
    bool PreservesOrientation() const 
    {
        return Determinant() > 0;
    }

    // Checking whether a transformation is isometric
    bool IsIsometry() const 
    {
        // For isometry: A^T * A = I
        double a = Matrix11 * Matrix11 + Matrix21 * Matrix21;
        double b = Matrix11 * Matrix12 + Matrix21 * Matrix22;
        double c = Matrix12 * Matrix12 + Matrix22 * Matrix22;

        return std::abs(a - 1.0) < 1e-10 &&
            std::abs(b) < 1e-10 &&
            std::abs(c - 1.0) < 1e-10;
    }

    // Getting Transformation Components
    void GetMatrix(double& a11, double& a12, double& a21, double& a22) const 
    {
        a11 = Matrix11; a12 = Matrix12; a21 = Matrix21; a22 = Matrix22;
    }

    void GetTranslation(double& tx, double& ty) const 
    {
        tx = Matrix13; ty = Matrix23;
    }

    // Decomposition of a transformation into components
    void Decompose(Complex& translation, Complex& scale, double& rotation_angle) const 
    {
        // Translation
        translation = Complex(Matrix13, Matrix23);

        // Scale and rotation from a 2x2 matrix
        double scale_x = std::sqrt(Matrix11 * Matrix11 + Matrix21 * Matrix21);
        double scale_y = std::sqrt(Matrix12 * Matrix12 + Matrix22 * Matrix22);
        scale = Complex(scale_x, scale_y);

        // Rotation angle
        rotation_angle = std::atan2(Matrix21, Matrix11);
    }

    // String representation
    std::string ToString() const 
    {
        return "AffineTransform([" +
            std::to_string(Matrix11) + ", " + std::to_string(Matrix12) + ", " + std::to_string(Matrix13) + "], [" +
            std::to_string(Matrix21) + ", " + std::to_string(Matrix22) + ", " + std::to_string(Matrix23) + "])";
    }

    // Comparison Operators
    bool operator==(const AffineTransform& other) const 
    {
        return std::abs(Matrix11 - other.Matrix11) < 1e-10 &&
            std::abs(Matrix12 - other.Matrix12) < 1e-10 &&
            std::abs(Matrix13 - other.Matrix13) < 1e-10 &&
            std::abs(Matrix21 - other.Matrix21) < 1e-10 &&
            std::abs(Matrix22 - other.Matrix22) < 1e-10 &&
            std::abs(Matrix23 - other.Matrix23) < 1e-10;
    }

    bool operator!=(const AffineTransform& other) const 
    {
        return !(*this == other);
    }
};

// Inference operator
std::ostream& operator<<(std::ostream& Output, const AffineTransform& Transform) 
{
    Output << Transform.ToString();
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

    // Testing Homothety class
    std::cout << "\n=== Testing Homothety ===" << "\n";

    // Create a homothety with center at (1,1) and coefficient 2
    Homothety H1(Complex(1.0, 1.0), 2.0);
    std::cout << "H1 = " << H1 << "\n";

    // Test applying homothety to a point
    Complex point(3.0, 2.0);
    Complex transformed = H1.ApplyTo(point);
    std::cout << "H1.applyTo(" << point << ") = " << transformed << "\n";

    // Test composition of homotheties
    Homothety H2(Complex(0.0, 0.0), 0.5);
    Homothety H3 = H1.ComposeWith(H2);
    std::cout << "H1.composeWith(H2) = " << H3 << "\n";

    // Test inverse homothety
    Homothety H1_inv = H1.getInverse();
    std::cout << "H1.getInverse() = " << H1_inv << "\n";

    // Test power of homothety
    Homothety H1_pow = H1.raisedTo(2);
    std::cout << "H1.raisedTo(2) = " << H1_pow << "\n";

    // Test creating homothety from point pairs
    Complex A(1.0, 1.0), A1(3.0, 3.0);  // A -> A1
    Complex B(2.0, 0.0), B1(4.0, -2.0); // B -> B1

    Homothety H4 = Homothety::fromTwoPairs(A, A1, B, B1);
    std::cout << "Homothety from pairs: " << H4 << "\n";

    // Verify the homothety works correctly
    std::cout << "H4.applyTo(A) = " << H4.ApplyTo(A) << " (should be " << A1 << ")" << "\n";
    std::cout << "H4.applyTo(B) = " << H4.ApplyTo(B) << " (should be " << B1 << ")" << "\n";

    // Testing affine transformations
    std::cout << "\n=== Testing Affine Transformations ===" << "\n";

    // Creating Various Transformations
    AffineTransform translate = AffineTransform::Translation(2.0, 3.0);
    AffineTransform scale = AffineTransform::Scaling(2.0, 1.5);
    AffineTransform rotate = AffineTransform::RotationDegrees(45.0);

    std::cout << "Translation: " << translate << "\n";
    std::cout << "Scaling: " << scale << "\n";
    std::cout << "Rotation: " << rotate << "\n";

    // Applying transformations to a point
    Complex point(1.0, 1.0);
    Complex translated_point = translate.ApplyTo(point);
    Complex scaled_point = scale.ApplyTo(point);
    Complex rotated_point = rotate.ApplyTo(point);

    std::cout << "Point: " << point << "\n";
    std::cout << "Translated: " << translated_point << "\n";
    std::cout << "Scaled: " << scaled_point << "\n";
    std::cout << "Rotated: " << rotated_point << "\n";

    // Composition of transformations
    AffineTransform composite = translate * scale * rotate;
    std::cout << "Composite transform: " << composite << "\n";

    Complex transformed_point = composite.ApplyTo(point);
    std::cout << "Composite applied: " << transformed_point << "\n";

    // Reverse conversion
    AffineTransform inverse = composite.Inverse();
    Complex original_point = inverse.ApplyTo(transformed_point);
    std::cout << "Inverse applied: " << original_point << " (should be " << point << ")\n";

    // Transformation around an arbitrary point
    Complex center(1.0, 1.0);
    AffineTransform rotate_around_center = AffineTransform::RotationDegrees(90.0, center);
    Complex rotated_around_center = rotate_around_center.ApplyTo(Complex(2.0, 1.0));
    std::cout << "Rotation around center " << center << ": " << rotated_around_center << "\n";

    // Transformation decomposition
    Complex trans, scl;
    double angle;
    composite.Decompose(trans, scl, angle);
    std::cout << "Decomposition - Translation: " << trans << ", Scale: " << scl << ", Angle: " << angle * 180.0 / PI << "°\n";

    return 0;
}
