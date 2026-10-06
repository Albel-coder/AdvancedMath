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
    double angleRadians;

public:

    OrientedAngle(double radians = 0.0) : angleRadians(radians) {}

    static OrientedAngle fromDegrees(double degrees);   
    static OrientedAngle fromRadians(double radians);

    OrientedAngle normalized() const;

    // Normalization to the range [-PI, PI)
    OrientedAngle normalizedSigned() const;    

    double radians() const;
    double degrees() const;

    // Arithmetic operations
    OrientedAngle operator+(const OrientedAngle& other) const;
    OrientedAngle operator-(const OrientedAngle& other) const;
    OrientedAngle operator*(double scalar) const;
    OrientedAngle operator/(double scalar) const;

    // Assignment Operators
    OrientedAngle& operator+=(const OrientedAngle& other);
    OrientedAngle& operator-=(const OrientedAngle& other);
    OrientedAngle& operator*=(double scalar);
    OrientedAngle& operator/=(double scalar);

    OrientedAngle operator+() const;
    OrientedAngle operator-() const;

    bool operator==(const OrientedAngle& other) const;
    bool operator!=(const OrientedAngle& other) const;
    bool operator<(const OrientedAngle& other) const;
    bool operator<=(const OrientedAngle& other) const;
    bool operator>(const OrientedAngle& other) const;
    bool operator>=(const OrientedAngle& other) const;

    double sin() const;
    double cos() const;
    double tan() const;

    static OrientedAngle arcsin(double value);
    static OrientedAngle arccos(double value);
    static OrientedAngle arctan(double value);
    static OrientedAngle arctan2(double y, double x);

    static OrientedAngle zero();
    static OrientedAngle right();
    static OrientedAngle straight();
    static OrientedAngle full();

    bool isZero() const;
    bool isRight() const;
    bool isStraight() const;
    bool isAcute() const;
    bool isObtuse() const;
    bool isReflex() const;

    // Additional angle (sum up to 90)
    OrientedAngle complementary() const;
    // Adjacent angle (sum up to 180)
    OrientedAngle supplementary() const;

    std::string toString() const;
    std::string toStringRadians() const;
};
