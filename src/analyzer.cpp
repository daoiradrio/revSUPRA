#include <analyzer.hpp>


Analyzer::Analyzer(){};


Analyzer::~Analyzer(){};


void Analyzer::remove_doubles(std::string filepath, std::string filename){
    int i, j, k, l;
    int counter;
    Structure struc1;
    Structure struc2;
    std::string file1;
    std::string file2;
    std::string command;
    std::vector<std::string> files;
    std::ifstream filestream;
    std::vector<std::vector<double>> cost_mat;
    std::vector<int> assignment;
    double element_term;
    double cost;
    Eigen::Vector3d diff_vec;
    Eigen::MatrixX3d matched_coords1;
    Eigen::MatrixX3d matched_coords2;

    if (filepath.back() != '/'){
        filepath = filepath + "/";
    }

    command = "ls " + filepath + " | grep " + filename + " > " + filepath + "files.tmp";
    system(command.c_str());
    filestream.open(filepath + "files.tmp");
    if (filestream.is_open()){
        while (getline(filestream, file1)){
            files.push_back(file1);
        }
    }
    else{
        std::cout << "ERROR IN ANALYZER MODULE: COULD NOT OPEN FILES.TMP" << std::endl;
    }
    filestream.close();
    command = "rm -f " + filepath + "files.tmp";
    system(command.c_str());

    counter = 0;
    struc1.read_xyz(filepath + files[0]);
    cost_mat.resize(struc1.n_atoms, std::vector<double>(struc1.n_atoms, 0.0));
    matched_coords1.resize(struc1.n_atoms, 3);
    matched_coords1.setZero();
    matched_coords2.resize(struc1.n_atoms, 3);
    matched_coords2.setZero();

    std::cout << "Removing duplicate structures..." << std::endl;

    for (i = 0; i < files.size()-1; i++){
        file1 = filepath + "/" + filename + std::to_string(i) + ".xyz";
        struc1.read_xyz(file1);
        for (j = i + 1; j < files.size(); j++){
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
                counter++;
                break;
            }
        }
    }

    std::cout << "Individual conformers: " << files.size()-counter << std::endl;

    return;
}


void Analyzer::extract_energies(std::string folderpath, std::string foldername){
    int i, j;
    int dummy;
    double energy;
    std::string command;
    std::string filepath;
    std::string line;
    std::vector<std::string> folders;
    std::ifstream filestream;

    if (folderpath.back() != '/'){
        folderpath = folderpath + "/";
    }

    command = "ls " + folderpath + " | grep " + foldername + " > " + folderpath + "folders.tmp";
    system(command.c_str());
    filestream.open(folderpath + "folders.tmp");
    if (filestream.is_open()){
        while (getline(filestream, line)){
            folders.push_back(line + "/");
        }
    }
    else{
        std::cout << "ERROR IN ANALYZER MODULE: COULD NOT OPEN FOLDERS.TMP" << std::endl;
    }
    filestream.close();
    command = "rm -f " + folderpath + "folders.tmp";
    system(command.c_str());

    for (i = 0; i < folders.size(); i++){
        filepath = folderpath + foldername + std::to_string(i) + "/uffenergy";
        filestream.open(filepath);
        filestream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        if (filestream.is_open()){
            getline(filestream, line);
            std::stringstream linestream(line);
            linestream >> dummy >> energy;
            this->container.push_back(std::make_pair(energy, i));
        }
        filestream.close();
	}

    //std::sort(this->container.begin(), this->container.end(), this->sort_func);
    std::sort(
        this->container.begin(), this->container.end(), 
        [](std::pair<double, int>& a, std::pair<double, int>& b){
            return (a.first < b.first);
        }
    );
}


void Analyzer::divide_and_conquer_remove_doubles(std::string filepath, std::string filename){
    if (this->container.size() == 0){
        std::cout << "\nNO ENERGIES FOUND, SWITCHING TO STANDARD DOUBLE STRUCTURE REMOVAL\n" << std::endl;
        return this->remove_doubles(filepath, filename);
    }

    int index1, index2;
    double energy, min_energy, max_energy;
    double element_term;
    double cost;
    std::string command;
    std::string file1, file2;
    std::ifstream filestream;
    std::vector<std::string> files;
    std::vector<double> item;
    std::vector<std::pair<double, int>> copy_container;
    std::vector<std::pair<double, int>>::iterator it;
    Structure struc1;
    Structure struc2;

    int counter = 0;
    double tol = 0.0001;

    if (filepath.back() != '/'){
        filepath = filepath + "/";
    }

    command = "ls " + filepath + " | grep " + filename + " > " + filepath + "files.tmp";
    system(command.c_str());
    filestream.open(filepath + "files.tmp");
    if (filestream.is_open()){
        while (getline(filestream, file1)){
            files.push_back(file1);
        }
    }
    else{
        std::cout << "ERROR IN ANALYZER MODULE: COULD NOT OPEN FILES.TMP" << std::endl;
    }
    filestream.close();
    command = "rm -f " + filepath + "files.tmp";
    system(command.c_str());

    for (std::pair<double, int> item: container){
        energy = item.first;
        index1 = item.second;
        file1 = filepath + files[index1];
        struc1.read_xyz(file1);
        copy_container = this->container;
        it = std::find_if(
            copy_container.begin(), copy_container.end(), 
            [&](std::pair<double, int>& item){
                return (item.first == energy);
            }
        );
        copy_container.erase(it);
        min_energy = energy - tol;
        max_energy = energy + tol;
        /*std::vector<std::pair<double, int>>::iterator l = std::lower_bound(
            copy_container.begin(), copy_container.end(), min_energy,
            [](const std::pair<double, int>& item, const double& ref){
                return (item.first < ref);
            }
        );
        std::vector<std::pair<double, int>>::iterator r = std::lower_bound(
            copy_container.begin(), copy_container.end(), max_energy,
            [](const std::pair<double, int>& item, const double& ref){
                return (item.first < ref);
            }
        );*/
        std::vector<std::pair<double, int>>::iterator l = std::find_if(
            copy_container.begin(), copy_container.end(),
            [&](const std::pair<double, int>& item){
                return (item.first >= min_energy);
            }
        );
        std::vector<std::pair<double, int>>::iterator r = std::find_if(
            copy_container.begin(), copy_container.end(),
            [&](const std::pair<double, int>& item){
                return (item.first >= max_energy);
            }
        );
        for (; l != r; l++){
            index2 = (*l).second;
            file2 = filepath + files[index2];
            struc2.read_xyz(file2);
            if (this->doubles(struc1, struc2)){
                counter++;
                break;
            }
        }
    }

    std::cout << "Individual conformers: " << files.size()-counter << std::endl;

    /*for (std::vector<double> item: this->container){
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
                    counter++;
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
    }*/

    //std::cout << "Von " << files.size() << " erzeugten Konformeren sind " << files.size()-counter << " individuelle Strukturen\n";
}


bool Analyzer::doubles(Structure& struc1, Structure& struc2){
    int i, j;
    int n_atoms = struc1.n_atoms;
    bool doubles = false;
    double element_term;
    double cost;
    std::vector<std::vector<double>> cost_mat;
    std::vector<int> assignment;
    Eigen::Vector3d diff_vec;
    Eigen::MatrixX3d matched_coords1;
    Eigen::MatrixX3d matched_coords2;

    cost_mat.resize(n_atoms, std::vector<double>(n_atoms, 0.0));
    matched_coords1.resize(n_atoms, 3);
    matched_coords1.setZero();
    matched_coords2.resize(n_atoms, 3);
    matched_coords2.setZero();

    for (i = 0; i < n_atoms; i++){
        for (j = 0; j < i+1; j++){
            if (struc1.atoms[i]->element == struc2.atoms[j]->element){
                element_term = 0.0;
            }
            else{
                element_term = 100.0;
            }
            diff_vec = struc1.coords.row(i) - struc1.coords.row(j);
            cost = diff_vec.dot(diff_vec) + element_term;
            cost_mat[i][j] = cost;
            cost_mat[j][i] = cost;
        }
    }
  	assignment = hungarian(cost_mat);
    for (i = 0; i < assignment.size(); i++){
        matched_coords1(i, 0) = struc1.coords(i, 0);
        matched_coords1(i, 1) = struc1.coords(i, 1);
        matched_coords1(i, 2) = struc1.coords(i, 2);
        matched_coords2(i, 0) = struc2.coords(assignment[i], 0);
        matched_coords2(i, 1) = struc2.coords(assignment[i], 1);
        matched_coords2(i, 2) = struc2.coords(assignment[i], 2);
    }
    if (this->rmsd(matched_coords1, matched_coords2) <= 0.1){
        doubles = true;
    }

    return doubles;
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
