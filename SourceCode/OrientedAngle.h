#include <stdexcept>
#include <string>
#include <vector>
#include <memory>
#include <tuple>
#include <cmath>
#include <algorithm>
#include <sstream>

const double PI = 3.1415926535897932;

class OrientedAngle
{
private:

    double angle_rad; // Angle in radians

public:

    // Constructors
    OrientedAngle(double rad = 0.0) : angle_rad(rad) {}

    // Creation from degrees
    static OrientedAngle fromDegrees(double degrees)
    {
        return OrientedAngle(degrees * PI / 180.0);
    }

    // Creation from radians
    static OrientedAngle fromRadians(double radians)
    {
        return OrientedAngle(radians);
    }

    // Normalize the angle to the range [0, 2PI)
    OrientedAngle normalized() const
    {
        double normalized = std::fmod(angle_rad, 2 * PI);
        if (normalized < 0)
        {
            normalized += 2 * PI;
        }
        return OrientedAngle(normalized);
    }

    // Normalization to the range [-PI, PI)
    OrientedAngle normalizedSigned() const
    {
        double normalized = std::fmod(angle_rad, 2 * PI);
        if (normalized > PI)
        {
            normalized -= 2 * PI;
        }
        else if (normalized <= -PI)
        {
            normalized += 2 * PI;
        }

        return OrientedAngle(normalized);
    }

    // Getting a value in radians
    double radians() const
    {
        return angle_rad;
    }

    // Getting the value in degrees
    double degrees() const
    {
        return angle_rad * 180.0 / PI;
    }

    // Arithmetic operations
    OrientedAngle operator+(const OrientedAngle& other) const
    {
        return OrientedAngle(angle_rad + other.angle_rad);
    }

    OrientedAngle operator-(const OrientedAngle& other) const
    {
        return OrientedAngle(angle_rad - other.angle_rad);
    }

    OrientedAngle operator*(double scalar) const
    {
        return OrientedAngle(angle_rad * scalar);
    }

    OrientedAngle operator/(double scalar) const
    {
        if (scalar == 0.0)
        {
            throw std::runtime_error("Division by zero");
        }
        return OrientedAngle(angle_rad / scalar);
    }

    // Assignment Operators
    OrientedAngle& operator+=(const OrientedAngle& other)
    {
        angle_rad += other.angle_rad;
        return *this;
    }

    OrientedAngle& operator-=(const OrientedAngle& other)
    {
        angle_rad -= other.angle_rad;
        return *this;
    }

    OrientedAngle& operator*=(double scalar)
    {
        angle_rad *= scalar;
        return *this;
    }

    OrientedAngle& operator/=(double scalar)
    {
        if (scalar == 0.0)
        {
            throw std::runtime_error("Division by zero");
        }
        angle_rad /= scalar;
        return *this;
    }

    // Unary operators
    OrientedAngle operator+() const
    {
        return *this;
    }

    OrientedAngle operator-() const
    {
        return OrientedAngle(-angle_rad);
    }

    // Comparison Operators
    bool operator==(const OrientedAngle& other) const
    {
        auto norm1 = this->normalized();
        auto norm2 = other.normalized();
        return std::abs(norm1.angle_rad - norm2.angle_rad) < 1e-10;
    }

    bool operator!=(const OrientedAngle& other) const
    {
        return !(*this == other);
    }

    bool operator<(const OrientedAngle& other) const
    {
        auto norm1 = this->normalized();
        auto norm2 = other.normalized();
        return norm1.angle_rad < norm2.angle_rad;
    }

    bool operator<=(const OrientedAngle& other) const
    {
        auto norm1 = this->normalized();
        auto norm2 = other.normalized();
        return norm1.angle_rad <= norm2.angle_rad;
    }

    bool operator>(const OrientedAngle& other) const
    {
        auto norm1 = this->normalized();
        auto norm2 = other.normalized();
        return norm1.angle_rad > norm2.angle_rad;
    }

    bool operator>=(const OrientedAngle& other) const
    {
        auto norm1 = this->normalized();
        auto norm2 = other.normalized();
        return norm1.angle_rad >= norm2.angle_rad;
    }

    // Trigonometric functions
    double sin() const
    {
        return std::sin(angle_rad);
    }

    double cos() const
    {
        return std::cos(angle_rad);
    }

    double tan() const
    {
        return std::tan(angle_rad);
    }

    // Inverse trigonometric functions (static methods)
    static OrientedAngle arcsin(double value)
    {
        if (value < -1.0 || value > 1.0)
        {
            throw std::runtime_error("Value out of range for arcsin");
        }
        return OrientedAngle(std::asin(value));
    }

    static OrientedAngle arccos(double value)
    {
        if (value < -1.0 || value > 1.0)
        {
            throw std::runtime_error("Value out of range for arccos");
        }

        return OrientedAngle(std::acos(value));
    }

    static OrientedAngle arctan(double value)
    {
        return OrientedAngle(std::atan(value));
    }

    static OrientedAngle arctan2(double y, double x)
    {
        return OrientedAngle(std::atan2(y, x));
    }

    // Special corners
    static OrientedAngle zero()
    {
        return OrientedAngle(0.0);
    }

    static OrientedAngle right()
    {
        return OrientedAngle(PI / 2.0);
    }

    static OrientedAngle straight()
    {
        return OrientedAngle(PI);
    }

    static OrientedAngle full()
    {
        return OrientedAngle(2 * PI);
    }

    // Checking special cases
    bool isZero() const
    {
        auto norm = this->normalized();
        return std::abs(norm.angle_rad) < 1e-10;
    }

    bool isRight() const
    {
        auto norm = this->normalized();
        return std::abs(norm.angle_rad - PI / 2.0) < 1e-10;
    }

    bool isStraight() const
    {
        auto norm = this->normalized();
        return std::abs(norm.angle_rad - PI) < 1e-10;
    }

    bool isAcute() const
    {
        auto norm = this->normalized();
        return norm.angle_rad > 0 && norm.angle_rad < PI / 2.0;
    }

    bool isObtuse() const
    {
        auto norm = this->normalized();
        return norm.angle_rad > PI / 2.0 && norm.angle_rad < PI;
    }

    bool isReflex() const
    {
        auto norm = this->normalized();
        return norm.angle_rad > PI && norm.angle_rad < 2 * PI;
    }

    // Additional angle (sum up to 90)
    OrientedAngle complementary() const
    {
        auto norm = this->normalized();
        if (norm.angle_rad > PI / 2.0)
        {
            throw std::runtime_error("Angle too large for complementary angle");
        }
        return OrientedAngle(PI / 2.0 - norm.angle_rad);
    }

    // Adjacent angle (sum up to 180)
    OrientedAngle supplementary() const
    {
        auto norm = this->normalized();
        if (norm.angle_rad > PI)
        {
            throw std::runtime_error("Angle too large for supplementary angle");
        }
        return OrientedAngle(PI - norm.angle_rad);
    }

    // String representation
    std::string toString() const
    {
        return std::to_string(degrees());
    }

    std::string toStringRadians() const
    {
        return std::to_string(angle_rad) + " rad";
    }
};

// Inference operator
std::ostream& operator<<(std::ostream& os, const OrientedAngle& angle)
{
    os << angle.toString();
    return os;
}

// Multiplication by a scalar (commutativity)
OrientedAngle operator*(double scalar, const OrientedAngle& angle)
{
    return angle * scalar;
}