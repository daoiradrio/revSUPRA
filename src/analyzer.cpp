#include <analyzer.hpp>


Analyzer::Analyzer(){};


Analyzer::~Analyzer(){};


void Analyzer::remove_doubles(std::string filepath, std::string filename, int n_files){
    int i, j;
    Structure struc1;
    Structure struc2;
    std::string file1;
    std::string file2;
    std::string command;

    for (i = 0; i < n_files-1; i++){
        file1 = filename + std::to_string(i) + ".xyz";
        struc1.read_xyz(file1);
        for (j = i + 1; j < n_files; j++){
            file2 = filename + std::to_string(j) + ".xyz";
            struc2.read_xyz(file2);
            //if (this->rmsd(struc1.coords, struc2.coords) <= 0.1){
            //    command = "rm " + filepath + "/" + file2;
            //    system(command.c_str());
            //}
        }
    }

    return;
}


double Analyzer::rmsd(Eigen::MatrixX3d coords1, Eigen::MatrixX3d coords2){
    int i;
    double det;
    double rmsd;
    Eigen::Matrix3Xd coords1_T;
    Eigen::Matrix3Xd coords2_T;
    Eigen::MatrixX3d H;
    Eigen::Matrix3d helper_mat;
    Eigen::Matrix3d R;
    Eigen::Vector3d center1;
    Eigen::Vector3d center2;
    
    center1 = {0.0, 0.0, 0.0};
    center2 = {0.0, 0.0, 0.0};
    for (i = 0; i < coords1.rows(); i++){
        center1 = center1 + coords1.row(i).transpose();
        center2 = center2 + coords2.row(i).transpose();
    }
    
    center1 = (1.0/(double)coords1.rows()) * center1;
    center2 = (1.0/(double)coords2.rows()) * center2;

    for (i = 0; i < coords1.rows(); i++){
        coords1.row(i) = coords1.row(i) - center1.transpose();
        coords2.row(i) = coords2.row(i) - center2.transpose();
    }
    
    coords1_T = coords1.transpose();

    H = coords1_T * coords2;

    Eigen::JacobiSVD<Eigen::MatrixX3d> svd(H, Eigen::ComputeFullV | Eigen::ComputeFullU);
    
    std::cout << svd.matrixU() << std::endl;
    std::cout << std::endl;
    std::cout << svd.matrixV() << std::endl;
    std::cout << std::endl;
    std::cout << svd.singularValues() << std::endl;
    std::cout << std::endl;

    det = (svd.matrixV(), svd.matrixU().transpose()).determinant();
    
    if (det >= 0.0){
        det = 1.0;
    }
    else{
        det = -1.0;
    }

    helper_mat << 1.0, 0.0, 0.0,
         0.0, 1.0, 0.0,
         0.0, 0.0, det;

    R = (svd.matrixV() * helper_mat) * svd.matrixU().transpose();

    std::cout << R << std::endl;
    std::cout << std::endl; 

    std::cout << coords2 << std::endl;
    std::cout << std::endl;

    for (i = 0; i < coords2.rows(); i++){
        coords2.row(i) = coords2.row(i) * R;
    }
    
    std::cout << coords2 << std::endl;
    std::cout << std::endl;

    rmsd = 0.0;
    for (i = 0; i < coords1.rows(); i++){
        rmsd = rmsd + (coords1.row(i) - coords2.row(i)).norm();
    }
    rmsd = (1.0/(double)coords1.rows()) * rmsd;

    std::cout << rmsd << std::endl;

    return rmsd;
}