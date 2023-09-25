//* Author: Zac Zhuo Zhang
//* Calculate the tangent circle center of two non-parallel lines in 2D-XY space (Z = 0).
//* The distance between the line-line intersection and the tangent point (d) is a given value.
//* Knowing that there are four tangent circles:
//* The outcome is determined by the sequence input points that defines the line directions.

#include <Eigen/Dense>

class TangentCircleCalculator2D
{
public:
    std::pair<Eigen::Vector3d, double> getTangentCircleCenter(const Eigen::Vector3d &l0p0, const Eigen::Vector3d &l0p1, const Eigen::Vector3d &l1p0, const Eigen::Vector3d &l1p1, double d, Eigen::MatrixXd &tanPoints);
    Eigen::Vector3d lineLineIntersection(const Eigen::Vector3d &l0p0, const Eigen::Vector3d &l0p1, const Eigen::Vector3d &l1p0, const Eigen::Vector3d &l1p1);
};

std::pair<Eigen::Vector3d, double> TangentCircleCalculator2D::getTangentCircleCenter(const Eigen::Vector3d &l0p0, const Eigen::Vector3d &l0p1, const Eigen::Vector3d &l1p0, const Eigen::Vector3d &l1p1, double d, Eigen::MatrixXd &tanPoints)
{
    //* Intersection
    Eigen::Vector3d intersection = lineLineIntersection(l0p0, l0p1, l1p0, l1p1);

    //* Calculate the radians angle between the two lines.
    double dotProduct = (l0p1 - l0p0).dot(l1p1 - l1p0);
    double length0 = (l0p1 - l0p0).norm();
    double length1 = (l1p1 - l1p0).norm();
    double angle = std::acos(dotProduct / (length0 * length1));

    //* Distance between the intersection point and the tangent point.
    double radius = d * std::tan(angle / 2);

    // * Get the center
    Eigen::Vector3d center = intersection - (d / std::cos(angle / 2)) * ((l0p1 - l0p0).normalized() + (l1p1 - l1p0).normalized()).normalized();

    tanPoints = Eigen::MatrixXd(2, 3);
    tanPoints.row(0) = intersection - (radius / (std::tan(angle / 2))) * ((l0p1 - l0p0).normalized());
    tanPoints.row(1) = intersection - (radius / (std::tan(angle / 2))) * ((l1p1 - l1p0).normalized());

    return std::pair<Eigen::Vector3d, double>(center, radius);
}

Eigen::Vector3d TangentCircleCalculator2D::lineLineIntersection(const Eigen::Vector3d &l0p0, const Eigen::Vector3d &l0p1, const Eigen::Vector3d &l1p0, const Eigen::Vector3d &l1p1)
{
    Eigen::Vector3d D0 = l0p1 - l0p0;
    Eigen::Vector3d D1 = l1p1 - l1p0;

    Eigen::Vector3d v = l1p0 - l0p0;
    Eigen::Vector3d crossD0D1 = D0.cross(D1);

    if (crossD0D1.norm() < 1e-6)
    {
        return Eigen::Vector3d::Zero();
    }

    double t = v.cross(D1).norm() / crossD0D1.norm();
    double s = v.cross(D0).norm() / crossD0D1.norm();

    Eigen::Vector3d intersectionPoint = l0p0 + t * D0;

    return intersectionPoint;
}
