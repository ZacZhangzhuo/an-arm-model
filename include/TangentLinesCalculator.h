//* Author: Zac Zhuo Zhang
//* Calculate the tangent line of two circles in 2D-XY space (Z = 0).
//* Knowing that there are four tangent lines for two non-intersecting circles:
//* getTangentLines0 calculates the two tangent lines that intersects.
//* getTangentLines1 calculates the other two tangent lines that do not intersect.

#include <Eigen/Dense>

class TangentLinesCalculator2D
{
public:
    Eigen::MatrixXd getTangentLines0(const Eigen::Vector3d &p0, double r0, const Eigen::Vector3d &p1, double r1);
    Eigen::MatrixXd getTangentLines1(const Eigen::Vector3d &p0, double r0, const Eigen::Vector3d &p1, double r1);
    Eigen::Vector3d Transition2D(double _x, double _y, double _angle, int fx);
};

Eigen::MatrixXd TangentLinesCalculator2D::getTangentLines0(const Eigen::Vector3d &p0, double r0, const Eigen::Vector3d &p1, double r1)
{
    Eigen::Vector3d r1LowPoint;
    Eigen::Vector3d r0UpPoint;
    Eigen::Vector3d r0LowPoint;
    Eigen::Vector3d r1UpPoint;

    double aa = std::atan2(p1.y() - p0.y(), p1.x() - p0.x());
    double centerLine = (p1 - p0).norm();
    double bb = std::acos((std::abs(r0 - r1)) / centerLine);
    double cc = bb - aa;


    double r2LowX = r1 * std::cos(cc);
    double r2LowY = r1 * std::sin(cc);
    double r1UpX = r0 * std::cos(cc);
    double r1UpY = r0 * std::sin(cc);

    r1LowPoint.x() = p1.x() - r2LowX;
    r1LowPoint.y() = p1.y() + r2LowY;
    r0UpPoint.x() = p0.x() + r1UpX;
    r0UpPoint.y() = p0.y() - r1UpY;

    Eigen::Vector3d r2Tempoint = Transition2D(r1LowPoint.x() - p1.x(), r1LowPoint.y() - p1.y(), 2 * bb, 1);
    r1UpPoint.x() = p1.x() + r2Tempoint.x();
    r1UpPoint.y() = p1.y() + r2Tempoint.y();

    Eigen::Vector3d r1Tempoint = Transition2D(r0UpPoint.x() - p0.x(), r0UpPoint.y() - p0.y(), 2 * bb, 1);
    r0LowPoint.x() = p0.x() + r1Tempoint.x();
    r0LowPoint.y() = p0.y() + r1Tempoint.y();

    Eigen::MatrixXd Points(4, 3);
    Points.row(0) = r0LowPoint;
    Points.row(1) = r1LowPoint;
    Points.row(2) = r0UpPoint;
    Points.row(3) = r1UpPoint;

    return Points;
}

Eigen::MatrixXd TangentLinesCalculator2D::getTangentLines1(const Eigen::Vector3d &p0, double r0, const Eigen::Vector3d &p1, double r1)
{
    Eigen::Vector3d r1LowPoint;
    Eigen::Vector3d r0UpPoint;
    Eigen::Vector3d r0LowPoint;
    Eigen::Vector3d r1UpPoint;


    double centerLine = (p1 - p0).norm();
    double radiusSum = r0 + r1;
    double aa = std::atan2(p1.y() - p0.y(), p1.x() - p0.x());
    double bb = std::acos(radiusSum / centerLine);

    double r0UpX = r0 * std::cos(aa + bb);
    double r0UpY = r0 * std::sin(aa + bb);
    r0UpPoint.x() = p0.x() + r0UpX;
    r0UpPoint.y() = p0.y() + r0UpY;

    double r1LowX = r1 * std::cos(aa - bb);
    double r1LowY = r1 * std::sin(aa - bb);
    r1LowPoint.x() = p1.x() + r1LowX;
    r1LowPoint.y() = p1.y() + r1LowY;


    double r1UpX = r1 * std::cos(aa + bb);
    double r1UpY = r1 * std::sin(aa + bb);
    r1UpPoint.x() = p1.x() + r1UpX;
    r1UpPoint.y() = p1.y() + r1UpY;

    double r0LowX = r0 * std::cos(aa - bb);
    double r0LowY = r0 * std::sin(aa - bb);
    r0LowPoint.x() = p0.x() - r0LowX;
    r0LowPoint.y() = p0.y() - r0LowY;

    Eigen::MatrixXd Points(4, 3);
    Points.row(0) = r0UpPoint;
    Points.row(1) = r1LowPoint;
    Points.row(2) = r1UpPoint;
    Points.row(3) = r0LowPoint;

    return Points;

}

Eigen::Vector3d TangentLinesCalculator2D::Transition2D(double _x, double _y, double _angle, int fx)
{
    double newX = 0;
    double newY = 0;

    if (fx == 1)
    {
        newX = _x * std::cos(_angle) - _y * std::sin(_angle);
        newY = _x * std::sin(_angle) + _y * std::cos(_angle);
    }
    else
    {
        newX = _x * std::cos(_angle) + _y * std::sin(_angle);
        newY = _y * std::cos(_angle) - _x * std::sin(_angle);
    }

    return Eigen::Vector3d(newX, newY, 0);
}
