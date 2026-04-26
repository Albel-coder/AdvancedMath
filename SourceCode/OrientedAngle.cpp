#include "OrientedAngle.h"

OrientedAngle OrientedAngle::fromDegrees(double degrees) {
    return OrientedAngle(degrees * PI / 180.0);
}

OrientedAngle OrientedAngle::fromRadians(double radians) {
    return OrientedAngle(radians);
}

OrientedAngle OrientedAngle::normalized() const {
    double normalized = std::fmod(angleRadians, 2 * PI);
    if (normalized < 0) {
        normalized += 2 * PI;
    }

    return OrientedAngle(normalized);
}

OrientedAngle OrientedAngle::normalizedSigned() const {
    double normalized = std::fmod(angleRadians, 2 * PI);
    if (normalized > PI) {
        normalized -= 2 * PI;
    }
    else if (normalized <= -PI) {
        normalized += 2 * PI;
    }

    return OrientedAngle(normalized);
}

double OrientedAngle::radians() const {
    return angleRadians;
}

double OrientedAngle::degrees() const {
    return angleRadians * 180.0 / PI;
}

OrientedAngle OrientedAngle::operator+(const OrientedAngle& other) const {
    return OrientedAngle(angleRadians + other.angleRadians);
}

OrientedAngle OrientedAngle::operator-(const OrientedAngle& other) const {
    return OrientedAngle(angleRadians - other.angleRadians);
}

OrientedAngle OrientedAngle::operator*(double scalar) const {
    return OrientedAngle(angleRadians * scalar);
}

OrientedAngle OrientedAngle::operator/(double scalar) const {
    if (scalar == 0.0) {
        throw std::runtime_error("Division by zero");
    }

    return OrientedAngle(angleRadians / scalar);
}

OrientedAngle& OrientedAngle::operator+=(const OrientedAngle& other) {
    angleRadians += other.angleRadians;
    return *this;
}

OrientedAngle& OrientedAngle::operator-=(const OrientedAngle& other) {
    angleRadians -= other.angleRadians;
    return *this;
}

OrientedAngle& OrientedAngle::operator*=(double scalar) {
    angleRadians *= scalar;
    return *this;
}

OrientedAngle& OrientedAngle::operator/=(double scalar) {
    if (scalar == 0.0) {
        throw std::runtime_error("Division by zero");
    }
    angleRadians /= scalar;
    return *this;
}

OrientedAngle OrientedAngle::operator+() const{
    return *this;
}

OrientedAngle OrientedAngle::operator-() const {
    return OrientedAngle(-angleRadians);
}

bool OrientedAngle::operator==(const OrientedAngle& other) const {
    auto firstNormalized = this->normalized();
    auto secondNormalized = other.normalized();

    return std::abs(firstNormalized.angleRadians - secondNormalized.angleRadians) < 1e-10;
}

bool OrientedAngle::operator!=(const OrientedAngle& other) const {
    return !(*this == other);
}

double OrientedAngle::sin() const {
    return std::sin(angleRadians);
}

double OrientedAngle::cos() const {
    return std::cos(angleRadians);
}

double OrientedAngle::tan() const {
    return std::tan(angleRadians);
}

bool OrientedAngle::operator<(const OrientedAngle& other) const {
    auto firstNormalized = this->normalized();
    auto secondNormalized = other.normalized();

    return firstNormalized.angleRadians < secondNormalized.angleRadians;
}

bool OrientedAngle::operator<=(const OrientedAngle& other) const {
    auto firstNormalized = this->normalized();
    auto secondNormalized = other.normalized();

    return firstNormalized.angleRadians <= secondNormalized.angleRadians;
}

bool OrientedAngle::operator>(const OrientedAngle& other) const {
    auto firstNormalized = this->normalized();
    auto secondNormalized = other.normalized();

    return firstNormalized.angleRadians > secondNormalized.angleRadians;
}

bool OrientedAngle::operator>=(const OrientedAngle& other) const {
    auto firstNormalized = this->normalized();
    auto secondNormalized = other.normalized();

    return firstNormalized.angleRadians >= secondNormalized.angleRadians;
}

OrientedAngle OrientedAngle::arcsin(double value) {
    if (value < -1.0 || value > 1.0) {
        throw std::runtime_error("Value out of range for arcsin");
    }

    return OrientedAngle(std::asin(value));
}

OrientedAngle OrientedAngle::arccos(double value) {
    if (value < -1.0 || value > 1.0) {
        throw std::runtime_error("Value out of range for arccos");
    }

    return OrientedAngle(std::acos(value));
}

OrientedAngle OrientedAngle::arctan(double value) {
    return OrientedAngle(std::atan(value));
}

OrientedAngle OrientedAngle::arctan2(double y, double x) {
    return OrientedAngle(std::atan2(x, y));
}

OrientedAngle OrientedAngle::zero() {
    return OrientedAngle(0.0);
}

OrientedAngle OrientedAngle::right() {
    return OrientedAngle(PI / 2.0);
}

OrientedAngle OrientedAngle::straight() {
    return OrientedAngle(PI);
}

OrientedAngle OrientedAngle::full() {
    return OrientedAngle(2 * PI);
}

bool OrientedAngle::isZero() const {
    auto normalized = this->normalized();
    return std::abs(normalized.angleRadians) < 1e-10;
}

bool OrientedAngle::isRight() const {
    auto normalized = this->normalized();
    return std::abs(normalized.angleRadians - PI / 2.0) < 1e-10;
}

bool OrientedAngle::isStraight() const {
    auto normalized = this->normalized();
    return std::abs(normalized.angleRadians - PI) < 1e-10;
}

bool OrientedAngle::isAcute() const {
    auto normalized = this->normalized();
    return normalized.angleRadians > 0 && normalized.angleRadians < PI / 2.0;
}

bool OrientedAngle::isObtuse() const {
    auto normalized = this->normalized();
    return normalized.angleRadians > PI / 2.0 && normalized.angleRadians < PI;
}

bool OrientedAngle::isReflex() const {
    auto normalized = this->normalized();
    return normalized.angleRadians > PI && normalized.angleRadians < 2 * PI;
}

OrientedAngle OrientedAngle::complementary() const {
    auto normalized = this->normalized();
    if (normalized.angleRadians > PI / 2.0) {
        throw std::runtime_error("Angle too large for complementary angle");
    }

    return OrientedAngle(PI / 2.0 - normalized.angleRadians);
}

OrientedAngle OrientedAngle::supplementary() const {
    auto normalized = this->normalized();
    if (normalized.angleRadians > PI) {
        throw std::runtime_error("Angle too large for supplementary angle");
    }

    return OrientedAngle(PI - normalized.angleRadians);
}

std::string OrientedAngle::toString() const {
    return std::to_string(degrees());
}

std::string OrientedAngle::toStringRadians() const {
    return std::to_string(angleRadians) + " radians";
}

std::ostream& operator<<(std::ostream& outputStream, const OrientedAngle& angle) {
    outputStream << angle.toString();
    return outputStream;
}

OrientedAngle operator*(double scalar, const OrientedAngle& angle) {
    return angle * scalar;
}
