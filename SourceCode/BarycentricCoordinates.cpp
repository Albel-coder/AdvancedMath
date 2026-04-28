#include "BarycentricCoordinates.h"

BarycentricCoordinates BarycentricCoordinates::fromCartesian(const Complex& point, const Complex& alphaValue, const Complex& betaValue, const Complex& gammaValue) {
    Complex value0 = betaValue - alphaValue;
    Complex value1 = gammaValue - alphaValue;
    Complex value2 = point - alphaValue;

    double d00 = value0.getReal() * value0.getReal() + value0.getImag() * value0.getImag();
    double d01 = value0.getReal() * value1.getReal() + value0.getImag() * value1.getImag();
    double d11 = value1.getReal() * value1.getReal() + value1.getImag() * value1.getImag();
    double d20 = value2.getReal() * value0.getReal() + value2.getImag() * value0.getImag();
    double d21 = value2.getReal() * value1.getReal() + value2.getImag() * value1.getImag();

    double denominator = d00 * d11 - d01 * d01;

    if (std::abs(denominator) < 1e-10) {
        throw std::runtime_error("Degenerate triangle in barycentric coordinates");
    }

    double betaResult = (d11 * d20 - d01 * d21) / denominator;
    double gammaResult = (d00 * d21 - d01 * d20) / denominator;
    double alphaResult = 1.0 - betaResult - gammaResult;

    return BarycentricCoordinates(alphaResult, betaResult, gammaResult);
}

BarycentricCoordinates BarycentricCoordinates::vertexA() {
    return BarycentricCoordinates(1.0, 0.0, 0.0);
}

BarycentricCoordinates BarycentricCoordinates::vertexB() {
    return BarycentricCoordinates(0.0, 1.0, 0.0);
}

BarycentricCoordinates BarycentricCoordinates::vertexC() {
    return BarycentricCoordinates(0.0, 0.0, 1.0);
}

BarycentricCoordinates BarycentricCoordinates::centroid() {
    return BarycentricCoordinates(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0);
}

BarycentricCoordinates BarycentricCoordinates::incenter() {
    // All coordinates are equal for an equilateral triangle
    return centroid();
}

double BarycentricCoordinates::getAlpha() const {
    return alpha;
}

double BarycentricCoordinates::getBeta() const {
    return beta;
}

double BarycentricCoordinates::getGamma() const {
    return gamma;
}

void BarycentricCoordinates::setAlpha(double alphaValue) {
    alpha = alphaValue;
}

void BarycentricCoordinates::setBeta(double betaValue) {
    beta = betaValue;
}

void BarycentricCoordinates::setGamma(double gammaValue) {
    gamma = gammaValue;
}

BarycentricCoordinates BarycentricCoordinates::normalized() const {
    double sum = alpha + beta + gamma;
    if (std::abs(sum) < 1e-10) {
        throw std::runtime_error("Cannot normalize barycentric coordinates with zero sum");
    }

    return BarycentricCoordinates(alpha / sum, beta / sum, gamma / sum);
}

bool BarycentricCoordinates::isValid() const {
    return alpha >= 0.0 && beta >= 0.0 && gamma >= 0.0;
}

bool BarycentricCoordinates::isInsideTriangle() const {
    return isValid() && std::abs(alpha + beta + gamma - 1.0) < 1e-10;
}

bool BarycentricCoordinates::isOnEdge() const {
    return isValid() && (std::abs(alpha) < 1e-10 ||
        std::abs(beta) < 1e-10 ||
        std::abs(gamma) < 1e-10);
}

bool BarycentricCoordinates::isOnVertex() const {
    return (std::abs(alpha - 1.0) < 1e-10 && std::abs(beta) < 1e-10 && std::abs(gamma) < 1e-10) ||
        (std::abs(beta - 1.0) < 1e-10 && std::abs(alpha) < 1e-10 && std::abs(gamma) < 1e-10) ||
        (std::abs(gamma - 1.0) < 1e-10 && std::abs(alpha) < 1e-10 && std::abs(beta) < 1e-10);
}

Complex BarycentricCoordinates::toCartesian(const Complex& alphaValue, const Complex& betaValue, const Complex& gammaValue) const {
    return alphaValue * alpha + betaValue * beta + gammaValue * gamma;
}

BarycentricCoordinates BarycentricCoordinates::operator+(const BarycentricCoordinates& other) const {
    return BarycentricCoordinates(alpha + other.alpha,
        beta + other.beta,
        gamma + other.gamma);
}

BarycentricCoordinates BarycentricCoordinates::operator-(const BarycentricCoordinates& other) const {
    return BarycentricCoordinates(alpha - other.alpha,
        beta - other.beta,
        gamma - other.gamma);
}

BarycentricCoordinates BarycentricCoordinates::operator*(double scalar) const {
    return BarycentricCoordinates(alpha * scalar,
        beta * scalar,
        gamma * scalar);
}

BarycentricCoordinates BarycentricCoordinates::operator/(double scalar) const  {
    if (std::abs(scalar) < 1e-10) {
        throw std::runtime_error("Division by zero in barycentric coordinates");
    }

    return BarycentricCoordinates(alpha / scalar, beta / scalar, gamma / scalar);
}

BarycentricCoordinates& BarycentricCoordinates::operator+=(const BarycentricCoordinates& other) {
    alpha += other.alpha;
    beta += other.beta;
    gamma += other.gamma;
    return *this;
}

BarycentricCoordinates& BarycentricCoordinates::operator-=(const BarycentricCoordinates& other) {
    alpha -= other.alpha;
    beta -= other.beta;
    gamma -= other.gamma;
    return *this;
}

BarycentricCoordinates& BarycentricCoordinates::operator*=(double scalar) {
    alpha *= scalar;
    beta *= scalar;
    gamma *= scalar;
    return *this;
}

BarycentricCoordinates& BarycentricCoordinates::operator/=(double scalar) {
    if (std::abs(scalar) < 1e-10) {
        throw std::runtime_error("Division by zero in barycentric coordinates");
    }
    alpha /= scalar;
    beta /= scalar;
    gamma /= scalar;
    return *this;
}

std::vector<double> BarycentricCoordinates::interpolateColor(const BarycentricCoordinates& coordinates, const std::vector<double>& colorAlpha, const std::vector<double>& colorBeta, const std::vector<double>& colorGamma) {
    if (colorAlpha.size() != colorBeta.size() || colorAlpha.size() != colorGamma.size()) {
        throw std::runtime_error("Color vectors must have same size");
    }

    std::vector<double> result(colorAlpha.size());
    for (size_t i = 0; i < colorAlpha.size(); ++i) {
        result[i] = coordinates.alpha * colorAlpha[i] + coordinates.beta * colorBeta[i] + coordinates.gamma * colorGamma[i];
    }

    return result;
}

double BarycentricCoordinates::distanceTo(const BarycentricCoordinates& other) const {
    double da = alpha - other.alpha;
    double db = beta - other.beta;
    double dg = gamma - other.gamma;

    return std::sqrt(da * da + db * db + dg * dg);
}

BarycentricCoordinates BarycentricCoordinates::lerp(const BarycentricCoordinates& alpha, const BarycentricCoordinates& beta, double theta) {
    return alpha * (1.0 - theta) + beta * theta;
}

bool BarycentricCoordinates::operator==(const BarycentricCoordinates& other) const {
    return std::abs(alpha - other.alpha) < 1e-10 &&
        std::abs(beta - other.beta) < 1e-10 &&
        std::abs(gamma - other.gamma) < 1e-10;
}

bool BarycentricCoordinates::operator!=(const BarycentricCoordinates& other) const {
    return !(*this == other);
}

std::string BarycentricCoordinates::toString() const {
    return "Barycentric(" +
        std::to_string(alpha) + ", " +
        std::to_string(beta) + ", " +
        std::to_string(gamma) + ")";
}

bool BarycentricCoordinates::isCentroid() const {
    auto normalized = this->normalized();
    return std::abs(normalized.alpha - 1.0 / 3.0) < 1e-10 &&
        std::abs(normalized.beta - 1.0 / 3.0) < 1e-10 &&
        std::abs(normalized.gamma - 1.0 / 3.0) < 1e-10;
}

std::ostream& operator<<(std::ostream& outputStream, const BarycentricCoordinates& barycentricCoordinates) {
    outputStream << barycentricCoordinates.toString();
    return outputStream;
}

template<typename T>
inline T BarycentricCoordinates::interpolate(const BarycentricCoordinates& coordinates, const T& alphaValue, const T& betaValue, const T& gammaValue) {
    return alphaValue * coordinates.alpha + betaValue * coordinates.beta + gammaValue * coordinates.gamma;
}
