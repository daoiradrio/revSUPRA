#include <iostream>
#include <vector>
#include <unordered_map>
#include <Eigen/Dense>


std::unordered_map<std::string, float> valence_radii({
    {"C", 0.75},
    {"N", 0.71},
    {"O", 0.66},
    {"H", 0.32},
    {"B", 0.85},
    {"F", 0.64},
    {"Cl", 0.99},
    {"Br", 1.14},
    {"I", 1.33}
});


std::unordered_map<std::string, int> max_valences({
    {"C", 4},
    {"N", 3},
    {"O", 2},
    {"B", 3},
    {"H", 1},
    {"F", 1},
    {"Cl", 1},
    {"Br", 1},
    {"I", 1},
});


typedef struct atom
{
    std::string label;
    std::vector<double> coords;
    std::vector<atom> connectivity;
}atom;


std::string get_element(std::string label);
//float kabsch_rmsd(std::vector<atom> atoms1, std::vector<atom> atoms2);


std::string get_element(std::string label)
{   
    std::string element;
    
    if (isalpha(label[1]))
    {
        element = std::string() + label[0] + label[1];
    }
    else
    {
        element = label[0];
    }

    return element;
}


/*
// source: GitHub Oleg Alexandrov (username oleg-alexandrov, project "eigen")
// https://github.com/oleg-alexandrov/projects/blob/master/eigen/Kabsch.cpp (26.03.2022)
// Description:
// "Given two sets of 3D points, find the rotation + translation + scale which best maps the first set to the second. 
// Source: http://en.wikipedia.org/wiki/Kabsch_algorithm"
Eigen::Affine3d Find3DAffineTransform(Eigen::Matrix3Xd matrix1, Eigen::Matrix3Xd matrix2)
//Eigen::Matrix3d Find3DAffineTransform(Eigen::Matrix3Xd matrix1, Eigen::Matrix3Xd matrix2)
{
    // Default output
    Eigen::Affine3d A;
    A.linear() = Eigen::Matrix3d::Identity(3, 3);
    A.translation() = Eigen::Vector3d::Zero();

    if (matrix1.cols() != matrix2.cols())
    {
        throw "Find3DAffineTransform(): input data mis-match";
    }

    // First find the scale, by finding the ratio of sums of some distances,
    // then bring the datasets to the same scale.
    double dist_matrix1 = 0;
    double dist_matrix2 = 0;
    for (int col = 0; col < matrix1.cols()-1; col++)
    {
        dist_matrix1  += (matrix1.col(col+1) - matrix1.col(col)).norm();
        dist_matrix2 += (matrix2.col(col+1) - matrix2.col(col)).norm();
    }
    //HIER MUSS ANGEPASST WERDEN UM RÜCKGABE-TYP EINZUHALTEN
    //if (dist_matrix1 <= 0 || dist_matrix2 <= 0)
    //{
    //    return A;
    //}
    double scale = dist_matrix2/dist_matrix1;
    matrix2 /= scale;

    // Find the centroids then shift to the origin
    Eigen::Vector3d matrix1_ctr = Eigen::Vector3d::Zero();
    Eigen::Vector3d matrix2_ctr = Eigen::Vector3d::Zero();
    for (int col = 0; col < matrix1.cols(); col++)
    {
        matrix1_ctr += matrix1.col(col);
        matrix2_ctr += matrix2.col(col);
    }
    matrix1_ctr /= matrix1.cols();
    matrix2_ctr /= matrix2.cols();
    for (int col = 0; col < matrix1.cols(); col++)
    {
        matrix1.col(col) -= matrix1_ctr;
        matrix2.col(col) -= matrix2_ctr;
    }

    // SVD
    Eigen::MatrixXd Cov = matrix1 * matrix2.transpose();
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Cov, Eigen::ComputeThinU | Eigen::ComputeThinV);

    // Find the rotation
    double d = (svd.matrixV() * svd.matrixU().transpose()).determinant();
    if (d > 0)
    {
        d = 1.0;
    }
    else
    {
        d = -1.0;
    }
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity(3, 3);
    I(2, 2) = d;
    Eigen::Matrix3d R = svd.matrixV() * I * svd.matrixU().transpose();

    // The final transform
    A.linear() = scale * R;
    A.translation() = scale*(matrix2_ctr - R*matrix1_ctr);
    
    //Eigen::Matrix3d rotation_matrix = A.rotation();
    //return rotation_matrix;
    return A;
}
*/


/*
float kabsch_rmsd(std::vector<atom> atoms1, std::vector<atom> atoms2)
{
    int n = atoms1.size();
    Eigen::MatrixXd coords1(n, 3);
    Eigen::MatrixXd coords2(n, 3);

    for (int i = 0; i < n; i++)
    {
        coords1.row(i) << atoms1[i].coords[0], atoms1[i].coords[1], atoms1[i].coords[2];
        coords2.row(i) << atoms2[i].coords[0], atoms2[i].coords[1], atoms2[i].coords[2];
    }

    Eigen::Vector3d coords1_ctr = Eigen::Vector3d::Zero();
    Eigen::Vector3d coords2_ctr = Eigen::Vector3d::Zero();
    for (int row = 0; row < coords1.rows(); row++)
    {
        coords1_ctr += coords1.row(row);
        coords2_ctr += coords2.row(row);
    }
    coords1_ctr /= coords1.rows();
    coords2_ctr /= coords2.rows();
    for (int row = 0; row < coords1.rows(); row++)
    {
        coords1.row(row) -= coords1_ctr;
        coords2.row(row) -= coords2_ctr;
    }

    Eigen::MatrixXd Cov = coords1.transpose() * coords2;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Cov, Eigen::ComputeFullU | Eigen::ComputeFullV);

    double d = (svd.matrixV().transpose() * svd.matrixU().transpose()).determinant();
    if (d > 0)
    {
        d = 1.0;
    }
    else
    {
        d = -1.0;
    }
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity(3, 3);
    I(2, 2) = d;

    Eigen::Matrix3d R = svd.matrixV().transpose() * I * svd.matrixU().transpose();

    coords2 = coords2 * R;

    float rmsd =  0.0;
    for (int row = 0; row < coords1.rows(); row++)
    {
        rmsd += pow(coords1.coeff(row, 0) - coords2.coeff(row, 0), 2) + 
                pow(coords1.coeff(row, 1) - coords2.coeff(row, 1), 2) + 
                pow(coords1.coeff(row, 2) - coords2.coeff(row, 2), 2);
    }
    rmsd = sqrt((1.0/(float(coords1.rows())) * rmsd));

    return rmsd;
}
*/


float quaternion_rmsd(std::vector<atom> atoms1, std::vector<atom> atoms2)
{
    int n = atoms1.size();
    Eigen::MatrixXd coords1(3, n);
    Eigen::MatrixXd coords2(3, n);

    for (int i = 0; i < n; i++)
    {
        coords1.col(i) << atoms1[i].coords[0], atoms1[i].coords[1], atoms1[i].coords[2];
        coords2.col(i) << atoms2[i].coords[0], atoms2[i].coords[1], atoms2[i].coords[2];
    }

    Eigen::Vector3d coords1_ctr = Eigen::Vector3d::Zero();
    Eigen::Vector3d coords2_ctr = Eigen::Vector3d::Zero();
    for (int col = 0; col < coords1.cols(); col++)
    {
        coords1_ctr += coords1.col(col);
        coords2_ctr += coords2.col(col);
    }
    coords1_ctr /= coords1.cols();
    coords2_ctr /= coords2.cols();
    for (int col = 0; col < coords1.cols(); col++)
    {
        coords1.col(col) -= coords1_ctr;
        coords2.col(col) -= coords2_ctr;
    }

    Eigen::MatrixXd R = coords1 * coords2.transpose();

    Eigen::Matrix4d F;
    F << R.coeff(0,0) + R.coeff(1,1) + R.coeff(2,2), R.coeff(1,2) - R.coeff(2,1), R.coeff(2,0) - R.coeff(0,2), R.coeff(0,1) - R.coeff(1,0),
         R.coeff(1,2) - R.coeff(2,1), R.coeff(0,0) - R.coeff(1,1) - R.coeff(2,2), R.coeff(0,1) + R.coeff(1,0), R.coeff(0,2) + R.coeff(2,0),
         R.coeff(2,0) - R.coeff(0,2), R.coeff(0,1) + R.coeff(1,0), -R.coeff(0,0) + R.coeff(1,1) - R.coeff(2,2), R.coeff(1,2) + R.coeff(2,1),
         R.coeff(0,1) - R.coeff(1,0), R.coeff(0,2) + R.coeff(2,0), R.coeff(1,2) + R.coeff(2,1), -R.coeff(0,0) - R.coeff(1,1) + R.coeff(2,2);
    
    Eigen::EigenSolver<Eigen::Matrix4d> eigensolver;
    eigensolver.compute(F);
    Eigen::VectorXd eigen_values = eigensolver.eigenvalues().real();
    Eigen::MatrixXd eigen_vectors = eigensolver.eigenvectors().real();

    double max = eigen_values[0];
    for (int i = 1; i < eigen_values.size(); i++)
    {
        if (max < eigen_values[i])
        {
            max = eigen_values[i];
        }
    }
    
    float rmsd = 0;
    for (int i = 0; i < coords1.cols(); i++)
    {
        rmsd += coords1.col(i).dot(coords1.col(i)) + coords2.col(i).dot(coords2.col(i));
    }
    rmsd -= 2 * max;
    rmsd = sqrt((1.0/(float)coords1.cols()) * rmsd);

    return rmsd;
}
