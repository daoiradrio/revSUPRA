#include <analyzer.hpp>


Analyzer::Analyzer(){};


Analyzer::~Analyzer(){};


void Analyzer::remove_doubles(std::string filepath, std::string filename, int n_files){
    std::cout << "Removing duplicate structures..." << std::endl;

    int i, j, k, l;
    Structure struc1;
    Structure struc2;
    std::string file1;
    std::string file2;
    std::string command;
    std::vector<std::vector<double>> cost_mat;
    std::vector<int> assignment;
    double element_term;
    double cost;
    Eigen::Vector3d diff_vec;
    Eigen::MatrixX3d matched_coords1;
    Eigen::MatrixX3d matched_coords2;

    struc1.read_xyz("../apply/" + filepath + "/" + filename + "0.xyz");
    cost_mat.resize(struc1.n_atoms, std::vector<double>(struc1.n_atoms, 0.0));
    matched_coords1.resize(struc1.n_atoms, 3);
    matched_coords1.setZero();
    matched_coords2.resize(struc1.n_atoms, 3);
    matched_coords2.setZero();

    for (i = 0; i < n_files-1; i++){
        file1 = filepath + "/" + filename + std::to_string(i) + ".xyz";
        struc1.read_xyz(file1);
        for (j = i + 1; j < n_files; j++){
            file2 = filepath + "/" + filename + std::to_string(j) + ".xyz";
            struc2.read_xyz(file2);
            for (k = 0; k < struc1.n_atoms; k++){
                for (l = 0; l < k+1; l++){
                    if (struc1.atoms[k]->element == struc2.atoms[l]->element){
                        element_term = 0.0;
                    }
                    else{
                        element_term = 100.0;
                    }
                    diff_vec = struc1.coords.row(k) - struc1.coords.row(l);
                    cost = diff_vec.dot(diff_vec) + element_term;
                    cost_mat[k][l] = cost;
                    cost_mat[l][k] = cost;
                }
            }
  	        assignment = hungarian(cost_mat);
            for (k = 0; k < assignment.size(); k++){
                matched_coords1(k, 0) = struc1.coords(k, 0);
                matched_coords1(k, 1) = struc1.coords(k, 1);
                matched_coords1(k, 2) = struc1.coords(k, 2);
                matched_coords2(k, 0) = struc2.coords(assignment[k], 0);
                matched_coords2(k, 1) = struc2.coords(assignment[k], 1);
                matched_coords2(k, 2) = struc2.coords(assignment[k], 2);
            }
            if (this->rmsd(matched_coords1, matched_coords2) <= 0.1){
                n_files -= 1;
                break;
            }
        }
    }

    std::cout << "Individual conformers: " << n_files << std::endl;

    return;
}


void Analyzer::extract_energies(std::string folderpath, std::string foldername, int n_folders){
    int i, j;
    int dummy;
    double energy;
    std::string filepath;
    std::string line;
    std::ifstream file;

    for (i = 0; i < n_folders; i++){
        filepath = folderpath + foldername + std::to_string(i) + "/uffenergy";
        file.open(filepath);
        file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        if (file.is_open()){
            getline(file, line);
            std::stringstream linestream(line);
            linestream >> dummy >> energy;
            this->container.push_back({energy, (double)i});
        }
        file.close();
	}

    std::sort(this->container.begin(), this->container.end(), this->sort_func);
}


void Analyzer::divide_and_conquer_remove_doubles(std::string filepath, std::string filename, int n_files){
    int i, j, k;
    int n;
    int l, m, r;
    int counter;
    double energy;
    double element_term;
    double cost;
    std::string file1, file2;
    std::vector<double> item;
    std::vector<std::vector<double>> copy_container;
    //std::vector<std::vector<double>>::iterator it;
    Structure struc1;
    Structure struc2;
    std::vector<std::vector<double>> cost_mat;
    std::vector<int> assignment;
    Eigen::Vector3d diff_vec;
    Eigen::MatrixX3d matched_coords1;
    Eigen::MatrixX3d matched_coords2;

    counter = n_files;
    struc1.read_xyz(filepath + filename + "0.xyz");
    cost_mat.resize(struc1.n_atoms, std::vector<double>(struc1.n_atoms, 0.0));
    matched_coords1.resize(struc1.n_atoms, 3);
    matched_coords1.setZero();
    matched_coords2.resize(struc1.n_atoms, 3);
    matched_coords2.setZero();

    for (std::vector<double> item: this->container){
        energy = item[0];
        n = (int)item[1];
        file1 = filepath + filename + std::to_string(n) + ".xyz";
        struc1.read_xyz(file1);
        copy_container = this->container;
        for (i = n_files-1; i >=0; i--){
            if (energy == copy_container[i][0]){
                copy_container.erase(copy_container.begin()+i);
            }
        }
        l = 0;
        r = copy_container.size()-1;
        while (l <= r){
            m = l + (r - l)/2;
            if (fabs(energy - copy_container[m][0]) < 0.0001){
                file2 = filepath + filename + std::to_string(m) + ".xyz";
                struc2.read_xyz(file2);
                for (j = 0; j < struc1.n_atoms; j++){
                    for (k = 0; k < j+1; k++){
                        if (struc1.atoms[j]->element == struc2.atoms[k]->element){
                            element_term = 0.0;
                        }
                        else{
                            element_term = 100.0;
                        }
                        diff_vec = struc1.coords.row(j) - struc1.coords.row(k);
                        cost = diff_vec.dot(diff_vec) + element_term;
                        cost_mat[j][k] = cost;
                        cost_mat[k][j] = cost;
                    }
                }
  	            assignment = hungarian(cost_mat);
                for (j = 0; j < assignment.size(); j++){
                    matched_coords1(j, 0) = struc1.coords(j, 0);
                    matched_coords1(j, 1) = struc1.coords(j, 1);
                    matched_coords1(j, 2) = struc1.coords(j, 2);
                    matched_coords2(j, 0) = struc2.coords(assignment[j], 0);
                    matched_coords2(j, 1) = struc2.coords(assignment[j], 1);
                    matched_coords2(j, 2) = struc2.coords(assignment[j], 2);
                }
                if (this->rmsd(matched_coords1, matched_coords2) <= 0.1){
                    counter -= 1;
                }
                break;
            }
            else if (energy < copy_container[m][0]){
                r = m - 1;
            }
            else{
                l = m + 1;
            }
        }
    }
    std::cout << "Von " << n_files << " erzeugten Konformeren sind " << counter << " individuelle Strukturen\n";
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
    
    /*std::cout << svd.matrixU() << std::endl;
    std::cout << std::endl;
    std::cout << svd.matrixV() << std::endl;
    std::cout << std::endl;
    std::cout << svd.singularValues() << std::endl;
    std::cout << std::endl;*/

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

    //std::cout << R << std::endl;
    //std::cout << std::endl; 

    //std::cout << coords2 << std::endl;
    //std::cout << std::endl;

    for (i = 0; i < coords2.rows(); i++){
        coords2.row(i) = coords2.row(i) * R;
    }
    
    //std::cout << coords2 << std::endl;
    //std::cout << std::endl;

    rmsd = 0.0;
    for (i = 0; i < coords1.rows(); i++){
        rmsd = rmsd + (coords1.row(i) - coords2.row(i)).norm();
    }
    rmsd = (1.0/(double)coords1.rows()) * rmsd;

    //std::cout << rmsd << std::endl;

    return rmsd;
}
