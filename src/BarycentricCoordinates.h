#include "Complex.h"

const double PI = 3.1415926535897932;

class BarycentricCoordinates
{
private:
    double alpha, beta, gamma;

public:

    BarycentricCoordinates(double alphaValue = 0.0, double betaValue = 0.0, double gammaValue = 0.0)
        : alpha(alphaValue), beta(betaValue), gamma(gammaValue) {
    }

    // Creation from Cartesian coordinates relative to a triangle
    static BarycentricCoordinates fromCartesian(const Complex& point, const Complex& alphaValue, const Complex& betaValue, const Complex& gammaValue);

    // Creating special points of a triangle
    static BarycentricCoordinates vertexA();
    static BarycentricCoordinates vertexB();
    static BarycentricCoordinates vertexC();

    static BarycentricCoordinates centroid();
    static BarycentricCoordinates incenter();

    double getAlpha() const;
    double getBeta() const;
    double getGamma() const;

    void setAlpha(double alphaValue);
    void setBeta(double betaValue);
    void setGamma(double gammaValue);

    // Normalization (sum of coordinates = 1)
    BarycentricCoordinates normalized() const;

    bool isValid() const;
    bool isInsideTriangle() const;
    bool isOnEdge() const;
    bool isOnVertex() const;
    // Transformation to Cartesian coordinates
    Complex toCartesian(const Complex& alphaValue, const Complex& betaValue, const Complex& gammaValue) const;

    BarycentricCoordinates operator+(const BarycentricCoordinates& other) const;
    BarycentricCoordinates operator-(const BarycentricCoordinates& other) const;
    BarycentricCoordinates operator*(double scalar) const;
    BarycentricCoordinates operator/(double scalar) const;
    BarycentricCoordinates& operator+=(const BarycentricCoordinates& other);
    BarycentricCoordinates& operator-=(const BarycentricCoordinates& other);
    BarycentricCoordinates& operator*=(double scalar);
    BarycentricCoordinates& operator/=(double scalar);

    template<typename T>
    static T interpolate(const BarycentricCoordinates& coordinates, const T& alphaValue, const T& betaValue, const T& gammaValue);

    static std::vector<double> interpolateColor(const BarycentricCoordinates& coordinates,
        const std::vector<double>& colorAlpha,
        const std::vector<double>& colorBeta,
        const std::vector<double>& colorGamma);

    double distanceTo(const BarycentricCoordinates& other) const;
    static BarycentricCoordinates lerp(const BarycentricCoordinates& alpha, const BarycentricCoordinates& beta, double theta);

    bool operator==(const BarycentricCoordinates& other) const;
    bool operator!=(const BarycentricCoordinates& other) const;

    std::string toString() const;
    bool isCentroid() const;
};

std::ostream& operator<<(std::ostream& outputStream, const BarycentricCoordinates& barycentricCoordinates);
