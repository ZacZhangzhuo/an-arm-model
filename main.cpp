//* Author: Zac Zhuo Zhang
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>

#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

#include "TangentLinesCalculator.h"
#include "TangentCircleCalculator.h"
#include "DivideArc.h"

Eigen::MatrixXd V;
Eigen::MatrixXi F;

TangentLinesCalculator2D tangentLinesCalculator;
TangentCircleCalculator2D tangentCircleCalculator;
CircleArcDivider2D circleArcDivider;

std::pair<Eigen::MatrixXd, Eigen::MatrixXi> ArmMeshGeneration(double R0, double R1, double R2, double Rad0, double Rad1, double length0, double length1, double d, double MeshMaxDistance, double thickness)
{

    //* GEOMETRIC MATH----------------------------------------------------------------------------------------------------------------------------------------------------------------
    Eigen::Vector3d P0(0, 0, 0);
    Eigen::Vector3d P1(length0 * std::cos(Rad0), length0 * std::sin(Rad0), 0);
    Eigen::Vector3d P2 = P1 + Eigen::Vector3d(std::cos(-Rad1) * length1, std::sin(-Rad1) * length1, 0);

    Eigen::MatrixXd tan0Points = tangentLinesCalculator.getTangentLines0(P0, R0, P1, R1);

    Eigen::MatrixXd tan1Points = tangentLinesCalculator.getTangentLines0(P1, R1, P2, R2);
    // std::cout << tan1Points.size() << std::endl;

    Eigen::MatrixXd arcPoints0 = circleArcDivider.divideArc2D(P0, R0, tan0Points.row(0), tan0Points.row(2), MeshMaxDistance);
    Eigen::MatrixXd arcPoints1 = circleArcDivider.divideArc2D(P1, R1, tan1Points.row(0), tan0Points.row(1), MeshMaxDistance);

    Eigen::MatrixXd tanPoints;
    std::pair<Eigen::Vector3d, double> centerRadius = tangentCircleCalculator.getTangentCircleCenter(tan0Points.row(2), tan0Points.row(3), tan1Points.row(3), tan1Points.row(2), d, tanPoints);
    Eigen::MatrixXd arcPoints2 = circleArcDivider.divideArc2D(centerRadius.first, centerRadius.second, tanPoints.row(1), tanPoints.row(0), MeshMaxDistance);

    Eigen::MatrixXd arcPoints3 = circleArcDivider.divideArc2D(P2, R2, tan1Points.row(3), tan1Points.row(1), MeshMaxDistance);


    // * CREATE MESH----------------------------------------------------------------------------------------------------------------------------------------------------------------
    // *Vertices Collection
    Eigen::MatrixXd P(2, 3);
    P.row(0) = P0;
    P.row(1) = P2;
    int indexP0 = arcPoints0.rows() + arcPoints1.rows() + arcPoints2.rows() + arcPoints3.rows();
    int indexP1 = arcPoints0.rows() + arcPoints1.rows() + arcPoints2.rows() + arcPoints3.rows() + 1;

    Eigen::MatrixXd frontVertices(arcPoints0.rows() + arcPoints1.rows() + arcPoints2.rows() + arcPoints3.rows() + P.rows(), 3);
    frontVertices << arcPoints0, arcPoints1, arcPoints2, arcPoints3, P;

    Eigen::MatrixXd backVertices(frontVertices.rows(), 3);
    for (int i = 0; i < backVertices.rows(); i++)
    {
        backVertices(i, 0) = frontVertices(i, 0);
        backVertices(i, 1) = frontVertices(i, 1);
        backVertices(i, 2) = -thickness;
    }
    Eigen::MatrixXd V(frontVertices.rows() + backVertices.rows(), 3);
    V << frontVertices, backVertices;

    // * Faces
    // Section 0 with its side
    Eigen::MatrixXi f0(arcPoints0.rows() - 1, 3);
    Eigen::MatrixXi f0s((arcPoints0.rows() - 1) * 2, 3);

    for (int i = 0; i < arcPoints0.rows() - 1; i++)
    {
        f0.row(i) = Eigen::Vector3i(i, i + 1, indexP0);
        f0s.row(i) = Eigen::Vector3i(i, frontVertices.rows() + i, frontVertices.rows() + i + 1);
        f0s.row(arcPoints0.rows() - 1 + i) = Eigen::Vector3i(i, frontVertices.rows() + i + 1, i + 1);
    }

    // Section 1 with its side
    Eigen::MatrixXi f1(3, 3);
    Eigen::MatrixXi f1s(4, 3);
    f1.row(0) = Eigen::Vector3i(arcPoints0.rows() - 1, arcPoints0.rows() + arcPoints1.rows() + arcPoints2.rows() - 1, indexP0);
    f1.row(1) = Eigen::Vector3i(indexP0, arcPoints0.rows() + arcPoints1.rows() + arcPoints2.rows() - 1, arcPoints0.rows() + arcPoints1.rows() - 1);
    f1.row(2) = Eigen::Vector3i(indexP0, arcPoints0.rows() + arcPoints1.rows() - 1, 0);

    f1s.row(0) = Eigen::Vector3i(arcPoints0.rows() - 1, frontVertices.rows() + arcPoints0.rows() + arcPoints1.rows() + arcPoints2.rows() - 1, arcPoints0.rows() + arcPoints1.rows() + arcPoints2.rows() - 1);
    f1s.row(1) = Eigen::Vector3i(arcPoints0.rows() - 1, frontVertices.rows() + arcPoints0.rows() - 1, frontVertices.rows() + arcPoints0.rows() + arcPoints1.rows() + arcPoints2.rows() - 1);
    f1s.row(2) = Eigen::Vector3i(0, frontVertices.rows() + arcPoints0.rows() + arcPoints1.rows() - 1, frontVertices.rows());
    f1s.row(3) = Eigen::Vector3i(0, arcPoints0.rows() + arcPoints1.rows() - 1, frontVertices.rows() + arcPoints0.rows() + arcPoints1.rows() - 1);

    // Section 2 with its side
    std::vector<Eigen::Vector3i> f2_vec;
    for (int i = 1; i < std::max(arcPoints1.rows(), arcPoints2.rows()); i++)
    {
        if (i < arcPoints1.rows())
            if (i < arcPoints2.rows())
                f2_vec.push_back(Eigen::Vector3i(arcPoints0.rows() + i - 1, arcPoints0.rows() + i, arcPoints0.rows() + arcPoints1.rows() + i - 1));
            else
                f2_vec.push_back(Eigen::Vector3i(arcPoints0.rows() + i - 1, arcPoints0.rows() + i, arcPoints0.rows() + arcPoints1.rows() + arcPoints2.rows() - 1));
        if (i < arcPoints2.rows())
            if (i < arcPoints1.rows())
                f2_vec.push_back(Eigen::Vector3i(arcPoints0.rows() + arcPoints1.rows() + i, arcPoints0.rows() + arcPoints1.rows() + i - 1, arcPoints0.rows() + i));
            else
                f2_vec.push_back(Eigen::Vector3i(arcPoints0.rows() + arcPoints1.rows() + i, arcPoints0.rows() + arcPoints1.rows() + i - 1, arcPoints0.rows() + arcPoints1.rows() - 1));
    }

    Eigen::MatrixXi f2(f2_vec.size(), 3);
    for (int i = 0; i < f2_vec.size(); ++i)
        f2.row(i) = f2_vec[i];

    Eigen::MatrixXi f2s0((arcPoints1.rows() - 1) * 2, 3);
    for (int i = 0; i < arcPoints1.rows() - 1; i++)
    {
        f2s0.row(i) = Eigen::Vector3i(arcPoints0.rows() + i, frontVertices.rows() + arcPoints0.rows() + i, frontVertices.rows() + arcPoints0.rows() + i + 1);
        f2s0.row(arcPoints1.rows() - 1 + i) = Eigen::Vector3i(arcPoints0.rows() + i, frontVertices.rows() + arcPoints0.rows() + i + 1, arcPoints0.rows() + i + 1);
    }
    Eigen::MatrixXi f2s1((arcPoints2.rows() - 1) * 2, 3);
    for (int i = 0; i < arcPoints2.rows() - 1; i++)
    {
        f2s1.row(i) = Eigen::Vector3i(arcPoints0.rows() + arcPoints1.rows() + i, frontVertices.rows() + arcPoints0.rows() + arcPoints1.rows() + i + 1, frontVertices.rows() + arcPoints0.rows() + arcPoints1.rows() + i);
        f2s1.row(arcPoints2.rows() - 1 + i) = Eigen::Vector3i(arcPoints0.rows() + arcPoints1.rows() + i, arcPoints0.rows() + arcPoints1.rows() + i + 1, frontVertices.rows() + arcPoints0.rows() + arcPoints1.rows() + i + 1);
    }

    // Section 3
    Eigen::MatrixXi f3(3, 3);
    Eigen::MatrixXi f3s(4, 3);

    f3.row(0) = Eigen::Vector3i(arcPoints0.rows() + arcPoints1.rows(), arcPoints0.rows() + arcPoints1.rows() + arcPoints2.rows(), indexP1);
    f3.row(1) = Eigen::Vector3i(arcPoints0.rows() + arcPoints1.rows(), indexP1, arcPoints0.rows());
    f3.row(2) = Eigen::Vector3i(arcPoints0.rows(), indexP1, arcPoints0.rows() + arcPoints1.rows() + arcPoints2.rows() + arcPoints3.rows() - 1);

    f3s.row(0) = Eigen::Vector3i(arcPoints0.rows() + arcPoints1.rows(), frontVertices.rows() + arcPoints0.rows() + arcPoints1.rows(), frontVertices.rows() + arcPoints0.rows() + arcPoints1.rows() + arcPoints2.rows());
    f3s.row(1) = Eigen::Vector3i(arcPoints0.rows() + arcPoints1.rows(), frontVertices.rows() + arcPoints0.rows() + arcPoints1.rows() + arcPoints2.rows(), arcPoints0.rows() + arcPoints1.rows() + arcPoints2.rows());
    f3s.row(2) = Eigen::Vector3i(arcPoints0.rows(), frontVertices.rows() + arcPoints0.rows() + arcPoints1.rows() + arcPoints2.rows() + arcPoints3.rows() - 1, frontVertices.rows() + arcPoints0.rows());
    f3s.row(3) = Eigen::Vector3i(arcPoints0.rows(), arcPoints0.rows() + arcPoints1.rows() + arcPoints2.rows() + arcPoints3.rows() - 1, frontVertices.rows() + arcPoints0.rows() + arcPoints1.rows() + arcPoints2.rows() + arcPoints3.rows() - 1);

    // Section4
    Eigen::MatrixXi f4(arcPoints3.rows() - 1, 3);
    Eigen::MatrixXi f4s((arcPoints3.rows() - 1) * 2, 3);
    for (int i = 0; i < arcPoints3.rows() - 1; i++)
    {
        f4.row(i) = Eigen::Vector3i(arcPoints0.rows() + arcPoints1.rows() + arcPoints2.rows() + i, arcPoints0.rows() + arcPoints1.rows() + arcPoints2.rows() + i + 1, indexP1);
        f4s.row(i) = Eigen::Vector3i(arcPoints0.rows() + arcPoints1.rows() + arcPoints2.rows() + i, frontVertices.rows() + arcPoints0.rows() + arcPoints1.rows() + arcPoints2.rows() + i, frontVertices.rows() + arcPoints0.rows() + arcPoints1.rows() + arcPoints2.rows() + i + 1);
        f4s.row(arcPoints3.rows() - 1 + i) = Eigen::Vector3i(arcPoints0.rows() + arcPoints1.rows() + arcPoints2.rows() + i, frontVertices.rows() + arcPoints0.rows() + arcPoints1.rows() + arcPoints2.rows() + i + 1, arcPoints0.rows() + arcPoints1.rows() + arcPoints2.rows() + i + 1);
    }

    // Back
    Eigen::MatrixXi frontFaces(f0.rows() + f1.rows() + f2.rows() + f3.rows() + f4.rows(), 3);
    frontFaces << f0, f1, f2, f3, f4;
    // std::cout << frontVertices.rows() << std::endl;
    // std::cout << V.rows() << std::endl;
    // std::cout << frontFaces.rows() << std::endl;
    Eigen::MatrixXi backFaces(frontFaces.rows(), 3);
    for (int i = 0; i < frontFaces.rows(); i++)
    {
        backFaces(i, 0) = frontFaces(i, 0) + frontVertices.rows();
        backFaces(i, 1) = frontFaces(i, 2) + frontVertices.rows();
        backFaces(i, 2) = frontFaces(i, 1) + frontVertices.rows();
    }
    // std::cout << backFaces << std::endl;

    // Concatenate
    Eigen::MatrixXi F(frontFaces.rows() + backFaces.rows() + f0s.rows() + f4s.rows() + f1s.rows() + f2s0.rows() + f2s1.rows() + f3s.rows(), 3);
    F << frontFaces, backFaces, f0s, f4s, f1s, f2s0, f2s1, f3s;

    return std::pair<Eigen::MatrixXd, Eigen::MatrixXi>(V, F);
}

int main(int argc, char *argv[])
{
    //* DEFAULT PARAMETERS (m) ----------------------------------------------------------------------------------------------------------------------------------------------------------------
    double R0 = 0.300;
    double R1 = 0.200;
    double R2 = 0.200;
    double Rad0 = M_PI / 6;
    double Rad1 = -M_PI / 6;
    double length0 = 1.200;
    double length1 = 1.700;
    double d = 0.200;
    double MeshMaxDistance = 0.050; // The resolution of the mesh
    double thickness = 0.50;

    //* DISPLAY----------------------------------------------------------------------------------------------------------------------------------------------------------------

    double Rmin = 0.1f;
    double Rmax = 0.4f;
    double Radmin = 0.1f;
    double Radmax = 1.50f;
    double lengthmin = 0.4f;
    double lengthmax = 3.0f;
    double dmin = 0.1f;
    double dmax = 0.5f;
    double MeshMaxDistancemin = 0.01f;
    double MeshMaxDistancemax = 0.1f;
    double thicknessmin = 0.1f;
    double thicknessmax = 1.0f;

    igl::opengl::glfw::Viewer viewer;

    std::pair<Eigen::MatrixXd, Eigen::MatrixXi> VF = ArmMeshGeneration(R0, R1, R2, Rad0, Rad1, length0, length1, d, MeshMaxDistance, thickness);
    viewer.data().set_mesh(VF.first, VF.second);

    igl::opengl::glfw::imgui::ImGuiPlugin plugin;
    viewer.plugins.push_back(&plugin);
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    plugin.widgets.push_back(&menu);

    menu.callback_draw_custom_window = [&]()
    {
        ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(250, 400), ImGuiCond_FirstUseEver);
        ImGui::Begin("The Arm", nullptr, ImGuiWindowFlags_NoSavedSettings);
        ImGui::PushItemWidth(-80);
        if (
            ImGui::DragScalar("R0", ImGuiDataType_Double, &R0, 0.1, &Rmin, &Rmax, "%.4f", (0.01F)) ||
            ImGui::DragScalar("R1", ImGuiDataType_Double, &R1, 0.1, &Rmin, &Rmax, "%.4f", (0.01F)) ||
            ImGui::DragScalar("R2", ImGuiDataType_Double, &R2, 0.1, &Rmin, &Rmax, "%.4f", (0.01F)) ||
            ImGui::DragScalar("Rad0", ImGuiDataType_Double, &Rad0, 0.1, &Radmin, &Radmax, "%.4f", (0.01F)) ||
            ImGui::DragScalar("Rad1", ImGuiDataType_Double, &Rad1, 0.1, &Radmin, &Radmax, "%.4f", (0.01F)) ||
            ImGui::DragScalar("length0", ImGuiDataType_Double, &length0, 0.1, &lengthmin, &lengthmax, "%.4f", (0.01F)) ||
            ImGui::DragScalar("length1", ImGuiDataType_Double, &length1, 0.1, &lengthmin, &lengthmax, "%.4f", (0.01F)) ||
            ImGui::DragScalar("d", ImGuiDataType_Double, &d, 0.1, &dmin, &dmax, "%.4f", (0.01F)) ||
            ImGui::DragScalar("MeshMaxDistance", ImGuiDataType_Double, &MeshMaxDistance, 0.1, &MeshMaxDistancemin, &MeshMaxDistancemax, "%.4f", (0.01F)) ||
            ImGui::DragScalar("thickness", ImGuiDataType_Double, &thickness, 0.1, &thicknessmin, &thicknessmax, "%.4f", (0.01F)))
        {
            viewer.data().clear();
            std::pair<Eigen::MatrixXd, Eigen::MatrixXi> VF = ArmMeshGeneration(R0, R1, R2, Rad0, Rad1, length0, length1, d, MeshMaxDistance, thickness);
            viewer.data().set_mesh(VF.first, VF.second);
        }

        ImGui::PopItemWidth();
        ImGui::End();
    };

    // *Launch
    viewer.launch();
}
