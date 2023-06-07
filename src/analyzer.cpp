#include <analyzer.hpp>



Analyzer::Analyzer(){};



Analyzer::~Analyzer(){};



void Analyzer::remove_doubles(std::string filepath, std::string filename){
    int i, j;
    int counter = 0;
    Structure struc1;
    Structure struc2;
    std::string file1;
    std::string file2;
    std::string command;
    std::vector<std::string> files;
    std::ifstream filestream;

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

    std::cout << "Removing duplicate structures..." << std::endl;

    int compares = 0;

    for (i = 0; i < files.size()-1; i++){
        file1 = filepath + files[i];
        struc1.read_xyz(file1);
        compares = 0;
        for (j = i + 1; j < files.size(); j++){
            compares++;
            file2 = filepath + files[j];
            struc2.read_xyz(file2);
            if (this->doubles(struc1, struc2)){
                counter++;
                break;
            }
        }
    }

    std::cout << "Individual conformers: " << files.size()-counter << std::endl;

    return;
}



/*void Analyzer::extract_energies(std::string folderpath, std::string foldername){
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
            this->energies.push_back(energy);
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

    return;
}*/



/*void Analyzer::divide_and_conquer_remove_doubles(std::string filepath, std::string filename){
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
    std::vector<std::pair<double, int>> working_container;
    std::vector<std::pair<double, int>>::iterator it1;
    std::vector<std::pair<double, int>>::iterator it2;
    Structure struc1;
    Structure struc2;

    int counter = 0;
    double tol = 0.025;

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

    int compares = 0;
    working_container = this->container;

    for (std::pair<double, int> item: container){
        energy = item.first;
        index1 = item.second;
        file1 = filepath + files[index1];
        struc1.read_xyz(file1);
        copy_container = working_container;
        //it = std::find_if(
        //    copy_container.begin(), copy_container.end(), 
        //    [&](std::pair<double, int>& item){
        //        return (item.first == energy);
        //    }
        //);
        for (it1 = copy_container.begin(), it2 = working_container.begin(); it1 != copy_container.end(); it1++, it2++){
            if ((*it1).first == energy){
                break;
            }
        }
        copy_container.erase(it1);
        min_energy = energy - tol;
        max_energy = energy + tol;
        std::vector<std::pair<double, int>>::iterator l = std::lower_bound(
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
        );
        //std::vector<std::pair<double, int>>::iterator l = std::find_if(
        //    copy_container.begin(), copy_container.end(),
        //    [&](const std::pair<double, int>& item){
        //        return (item.first >= min_energy);
        //    }
        //);
        //std::vector<std::pair<double, int>>::iterator r = std::find_if(
        //    copy_container.begin(), copy_container.end(),
        //    [&](const std::pair<double, int>& item){
        //        return (item.first >= max_energy);
        //    }
        //);
        compares = 0;
        for (; l != r; l++){
            compares++;
            index2 = (*l).second;
            file2 = filepath + files[index2];
            struc2.read_xyz(file2);
            if (this->doubles(struc1, struc2)){
                std::cout << index1 << " " << index2 << std::endl;
                std::cout << this->container[index1].first << " " << this->container[index2].first << std::endl;
                std::cout << compares << std::endl;
                std::cout << std::endl;
                working_container.erase(it2);
                counter++;
                break;
            }
        }
    }

    std::cout << "Individual conformers: " << files.size()-counter << std::endl;

    return;
}*/



bool Analyzer::doubles(Structure& struc1, Structure& struc2){
    int i, j;
    int n_atoms = struc1.n_atoms;
    bool doubles = false;
    double element_term;
    double cost;
    std::vector<std::vector<double>> cost_mat;
    std::vector<int> assignment;
    Eigen::Vector3d diff_vec;
    Eigen::MatrixX3d coords1 = struc1.coords;
    Eigen::MatrixX3d coords2 = struc2.coords;
    Eigen::MatrixX3d matched_coords1;
    Eigen::MatrixX3d matched_coords2;

    cost_mat.resize(n_atoms, std::vector<double>(n_atoms, 0.0));
    matched_coords1.resize(n_atoms, 3);
    matched_coords1.setZero();
    matched_coords2.resize(n_atoms, 3);
    matched_coords2.setZero();

    /*std::cout << "coords1" << std::endl;
    std::cout << struc1.coords << std::endl;
    std::cout << std::endl;*/

    this->kabsch(coords1, coords2);

    /*std::cout << "coords1" << std::endl;
    std::cout << struc1.coords << std::endl;
    std::cout << std::endl;*/

    for (i = 0; i < n_atoms; i++){
        for (j = 0; j < i+1; j++){
            if (struc1.atoms[i]->element == struc2.atoms[j]->element){
                element_term = 0.0;
            }
            else{
                element_term = 1000.0;
            }
            diff_vec = coords1.row(i) - coords2.row(j);
            cost = diff_vec.dot(diff_vec) + element_term;
            cost_mat[i][j] = cost;
            cost_mat[j][i] = cost;
        }
    }
  	assignment = hungarian(cost_mat);
    for (i = 0; i < assignment.size(); i++){
        matched_coords1(i, 0) = coords1(i, 0);
        matched_coords1(i, 1) = coords1(i, 1);
        matched_coords1(i, 2) = coords1(i, 2);
        matched_coords2(i, 0) = coords2(assignment[i], 0);
        matched_coords2(i, 1) = coords2(assignment[i], 1);
        matched_coords2(i, 2) = coords2(assignment[i], 2);
    }

    //if (this->rmsd(coords1, coords2) <= 0.1){
    if (this->rmsd(matched_coords1, matched_coords2) <= 0.1){
        doubles = true;
    }

    return doubles;
}



void Analyzer::kabsch(Eigen::MatrixX3d& coords1, Eigen::MatrixX3d& coords2)
{
    int i;
    double det;
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

    /*std::cout << "H:" << std::endl;
    std::cout << H << std::endl;
    std::cout << std::endl;*/

    //Eigen::JacobiSVD<Eigen::MatrixX3d> svd(H, Eigen::ComputeFullV | Eigen::ComputeFullU);
    Eigen::BDCSVD<Eigen::MatrixX3d> svd(H, Eigen::ComputeFullV | Eigen::ComputeFullU);
    
    /*std::cout << "U:" << std::endl;
    std::cout << svd.matrixU() << std::endl;
    std::cout << std::endl;

    std::cout << "S:" << std::endl;
    std::cout << svd.singularValues() << std::endl;
    std::cout << std::endl;

    std::cout << "V:" << std::endl;
    std::cout << svd.matrixV().transpose() << std::endl;
    std::cout << std::endl;*/

    det = (svd.matrixV() * svd.matrixU().transpose()).determinant();

    /*std::cout << "det:" << std::endl;
    std::cout << det << std::endl;
    std::cout << std::endl;*/
    
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

    /*std::cout << "R:" << std::endl;
    std::cout << R << std::endl;
    std::cout << std::endl;*/

    for (i = 0; i < coords2.rows(); i++){
        coords2.row(i) = coords2.row(i) * R;
    }

    return;
}



double Analyzer::rmsd(Eigen::MatrixX3d coords1, Eigen::MatrixX3d coords2)
{
    int             i;
    int             n = coords1.rows();
    double          rmsd = 0.0;
    Eigen::Vector3d diff_vec;

    /*std::cout << "coords1" << std::endl;
    std::cout << coords1 << std::endl;
    std::cout << std::endl;

    std::cout << "coords2 hinterher" << std::endl;
    std::cout << coords2 << std::endl;
    std::cout << std::endl;*/

    for (i = 0; i < n; i++){
        diff_vec = coords1.row(i) - coords2.row(i);
        rmsd    += diff_vec.dot(diff_vec);
    }
    rmsd = (1.0/(double)n) * rmsd;
    rmsd = sqrt(rmsd);

    //std::cout << "RMSD:" << std::endl;
    //std::cout << rmsd << std::endl;

    return rmsd;
}
