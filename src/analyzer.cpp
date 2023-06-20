#include <analyzer.hpp>



Analyzer::Analyzer(){};



Analyzer::~Analyzer(){};



void Analyzer::remove_doubles(
    std::string filepath, std::string filename, bool ignore_methyl, double rmsd_threshold
){
    int                         i, j;
    int                         counter = 0;
    std::string                 file1;
    std::string                 file2;
    std::string                 command;
    std::vector<std::string>    files;
    std::ifstream               filestream;
    Structure                   mol1;
    Structure                   mol2;

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

    for (i = 0; i < files.size()-1; i++){
        file1 = filepath + files[i];
        if (ignore_methyl){
            mol1.get_structure(file1);
        }
        else{
            mol1.read_xyz(file1);
        }
        for (j = i + 1; j < files.size(); j++){
            file2 = filepath + files[j];
            if (ignore_methyl){
                mol2.get_structure(file2);
            }
            else{
                mol2.read_xyz(file2);
            }
            if (this->doubles(file1, file2)){
                command = "rm " + file1;
                system(command.c_str());
                counter++;
                break;
            //if (this->doubles(mol1, mol2, rmsd_threshold, ignore_methyl)){
                /*if (mol1.energy && mol2.energy){
                    if (mol1.energy < mol2.energy){
                        std::cout << "Lösche mol2" << std::endl;
                        command = "rm " + file2;
                        system(command.c_str());
                        counter++;
                    }
                    else{
                        std::cout << "Lösche mol1" << std::endl;
                        command = "rm " + file1;
                        system(command.c_str());
                        counter++;
                        break;
                    }
                }
                else{
                    std::cout << "hier" << std::endl;
                    command = "rm " + file1;
                    system(command.c_str());
                    counter++;
                    break;
                }*/
            }
        }
    }

    std::cout << "Individual conformers: " << files.size()-counter << std::endl;

    return;
}



bool Analyzer::doubles(
    const Structure& mol1, const Structure& mol2, bool ignore_methyl, double rmsd_threshold
){
    int                                 i, j;
    bool                                doubles = false;
    Eigen::MatrixX3d                    coords1;
    Eigen::MatrixX3d                    coords2;
    std::vector<std::shared_ptr<Atom>>  atoms1;
    std::vector<std::shared_ptr<Atom>>  atoms2;

    if (ignore_methyl){
        this->remove_methly_atoms(mol1, coords1, atoms1);
        this->remove_methly_atoms(mol2, coords2, atoms2);
    }
    else{
        coords1 = mol1.coords;
        coords2 = mol2.coords;
        atoms1 = mol1.atoms;
        atoms2 = mol2.atoms;
    }

    this->kabsch(coords1, coords2);

    this->match_coords(atoms1, atoms2, coords1, coords2);

    if (this->rmsd(coords1, coords2) <= rmsd_threshold){
        doubles = true;
    }

    return doubles;
}



bool Analyzer::doubles(std::string file1, std::string file2, bool ignore_methyl, double rmsd_threshold)
{
    int                                 i, j;
    bool                                doubles = false;
    Structure                           mol1;
    Structure                           mol2;
    Eigen::MatrixX3d                    coords1;
    Eigen::MatrixX3d                    coords2;
    std::vector<std::shared_ptr<Atom>>  atoms1;
    std::vector<std::shared_ptr<Atom>>  atoms2;

    if (ignore_methyl){
        mol1.get_structure(file1);
        mol2.get_structure(file2);
        this->remove_methly_atoms(mol1, coords1, atoms1);
        this->remove_methly_atoms(mol2, coords2, atoms2);
    }
    else{
        mol1.read_xyz(file1);
        mol2.read_xyz(file2);
        coords1 = mol1.coords;
        coords2 = mol2.coords;
        atoms1 = mol1.atoms;
        atoms2 = mol2.atoms;
    }

    this->kabsch(coords1, coords2);

    this->match_coords(atoms1, atoms2, coords1, coords2);

    if (this->rmsd(coords1, coords2) <= rmsd_threshold){
        doubles = true;
    }

    return doubles;
}



void Analyzer::remove_methly_atoms(
        const Structure& mol,
        Eigen::MatrixX3d& coords,
        std::vector<std::shared_ptr<Atom>>& atoms
    )
{   
    int                             i, j;
    bool                            methyl_flag = false;
    std::vector<Eigen::Vector3d>    tmp_coords;

    for (i = 0; i < mol.atoms.size(); i++){
        methyl_flag = false;
        if (mol.atoms[i]->element == "C"){
            if (mol.atoms[i]->core_of_terminal_group){
                continue;
            }
        }
        else if (is_terminal_atom(mol.atoms[i]->element)){
            for (int j: mol.atoms[i]->bond_partners){
                if (mol.atoms[j]->element == "C"){
                    if (mol.atoms[j]->core_of_terminal_group){
                        methyl_flag = true;
                        break;
                    }
                }
            }
            if (methyl_flag){
                continue;
            }
        }
        atoms.push_back(mol.atoms[i]);
        tmp_coords.push_back(mol.coords.row(i));
    }
    coords.resize(tmp_coords.size(), 3);
    for (i = 0; i < tmp_coords.size(); i++){
        coords.row(i) = tmp_coords[i];
    }

    return;
}



void Analyzer::match_coords(
    std::vector<std::shared_ptr<Atom>> atoms1,
    std::vector<std::shared_ptr<Atom>> atoms2,
    Eigen::MatrixX3d coords1,
    Eigen::MatrixX3d& coords2
    )
{   
    int                                 i, j;
    int                                 n_atoms = atoms1.size();
    double                              element_term;
    double                              cost;
    std::vector<std::vector<double>>    cost_mat;
    std::vector<int>                    assignment;
    Eigen::Vector3d                     diff_vec;
    Eigen::MatrixX3d                    copy_coords2 = coords2;

    cost_mat.resize(n_atoms, std::vector<double>(n_atoms, 0.0));

    for (i = 0; i < n_atoms; i++){
        for (j = 0; j < i+1; j++){
            if (atoms1[i]->element == atoms2[j]->element){
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
        coords2(i, 0) = copy_coords2(assignment[i], 0);
        coords2(i, 1) = copy_coords2(assignment[i], 1);
        coords2(i, 2) = copy_coords2(assignment[i], 2);
    }

    return;
}



void Analyzer::kabsch(Eigen::MatrixX3d& coords1, Eigen::MatrixX3d& coords2)
{
    int                 i;
    double              det;
    Eigen::Matrix3Xd    coords1_T;
    Eigen::Matrix3Xd    coords2_T;
    Eigen::MatrixX3d    H;
    Eigen::Matrix3d     helper_mat;
    Eigen::Matrix3d     R;
    Eigen::Vector3d     center1;
    Eigen::Vector3d     center2;
    
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

    //Eigen::JacobiSVD<Eigen::MatrixX3d> svd(H, Eigen::ComputeFullV | Eigen::ComputeFullU);
    Eigen::BDCSVD<Eigen::MatrixX3d> svd(H, Eigen::ComputeFullV | Eigen::ComputeFullU);

    det = (svd.matrixV() * svd.matrixU().transpose()).determinant();
    
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

    for (i = 0; i < n; i++){
        diff_vec = coords1.row(i) - coords2.row(i);
        rmsd    += diff_vec.dot(diff_vec);
    }
    rmsd = (1.0/(double)n) * rmsd;
    rmsd = sqrt(rmsd);

    return rmsd;
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
