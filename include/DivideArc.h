//* Author: Zac Zhuo Zhang
//* Divide the Arc based on the given distance,
//* Shorter distance means smoother arc edges,
//* Longer distance means less points on the arc, minimum 2 points.

#include <Eigen/Dense>

class CircleArcDivider2D
{
public:
    Eigen::MatrixXd divideArc2D(const Eigen::Vector3d &center, double radius, const Eigen::Vector3d &pointA, const Eigen::Vector3d &pointB, double distance);
};

Eigen::MatrixXd CircleArcDivider2D::divideArc2D(const Eigen::Vector3d &center, double radius, const Eigen::Vector3d &pointA, const Eigen::Vector3d &pointB, double distance)
{

    //* Calculate the radians angle between the two points.
    double dotProduct = (pointA - center).dot(pointB - center);
    double lengthA = (pointA - center).norm();
    double lengthB = (pointB - center).norm();
    double rRange = std::acos(dotProduct / (lengthA * lengthB));
    // std::cout << "dotProduct / (lengthA * lengthB): " << dotProduct / (lengthA * lengthB) << std::endl;
    if (std::abs(dotProduct / (lengthA * lengthB)+1) < 0.000001)
        rRange = M_PI;

    // std::cout << "Radius: " << rRange << std::endl;

    //* Calculate the delta radius.
    double deltaRadius = 2 * std::asin(distance / (2 * radius));
    // std::cout << "Delta radius: " << deltaRadius << std::endl;

    //* Calculate the number of points.
    int numberOfPoints = std::ceil(rRange / deltaRadius) + 1;
    if (numberOfPoints < 2)
        numberOfPoints = 2;
    // std::cout << "Number of points: " << numberOfPoints << std::endl;
    Eigen::Vector3d vecA = pointA - center;
    Eigen::MatrixXd points(numberOfPoints, 3);

    points.row(0) = pointA;

    for (int i = 1; i < numberOfPoints - 1; i++)
    {
        double r = i * deltaRadius;

        Eigen::Matrix3d rotationMatrix;
        rotationMatrix << std::cos(r), -std::sin(r), 0,
            std::sin(r), std::cos(r), 0,
            0, 0, 1;

        Eigen::Vector3d rotatedVector = rotationMatrix * vecA;
        points.row(i) = rotatedVector + center;
    }

    points.row(numberOfPoints - 1) = pointB;

    // std::cout << "Number of points: " << points.size() << std::endl;

    return points;
}