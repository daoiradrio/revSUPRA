#include "conformergenerator.hpp"



ConformerGenerator::ConformerGenerator(Structure input_mol){
    this->mol = std::make_shared<Structure>(input_mol);
    this->workdir_name = "opt_dir";
    this->struc_filename = "struc.xyz";
    this->opt_struc_filename = "opt_struc.xyz";
};


ConformerGenerator::~ConformerGenerator(){};


void ConformerGenerator::generate_conformers(){
    this->get_torsions();

    /*
    // CHECK NUMBER OF CENTRAL TORSIONS BEFORE CYCLE DETECTION
    int atom1, atom2;
    for (bond torsion: this->central_torsions){
        atom1 = torsion.atom_index1;
        atom2 = torsion.atom_index2;
        std::cout << this->mol->atoms[atom1]->element << this->mol->atoms[atom1]->index << " " 
                  << this->mol->atoms[atom2]->element << this->mol->atoms[atom2]->index << std::endl;
    }
    */

    this->find_cycles();

    /*
    // CHECK NUMBER OF CENTRAL TORSIONS AFTER CYCLE DETECTION
    std::cout << "__________\n" << std::endl;
    for (bond torsion: this->central_torsions){
        atom1 = torsion.atom_index1;
        atom2 = torsion.atom_index2;
        std::cout << this->mol->atoms[atom1]->element << this->mol->atoms[atom1]->index << " " 
                  << this->mol->atoms[atom2]->element << this->mol->atoms[atom2]->index << std::endl;
    }
    */

    this->find_peptidebonds();

    // - ADJUST THE WHOLE THING FROM HERE SO THAT SELECTION RESTART REALLY STARTS AT THE VERY BEGINNING
    //   (IF PEPTITEBONDS SHOULD BE CONSIDERED)
    // - PACK INTO SEPARATE FUNCTION CALLED SELECTION MENU OR SO
    // *** START OF FUNCTION SELECTION MENU ***
    bool valid_mode_input = false;
    int mode_input;
    while (!valid_mode_input){
        std::cout << "Consider rotatable all bonds to terminal groups like -CH3, -NH2, -OH (1) "
                     "or ignore them (2) "
                     "or ignore just ignore bonds to methyl-like groups (3)?: ";
        std::cin >> mode_input;
        switch (mode_input){
            case 1:
                for (bond torsion: this->central_torsions){
                    this->torsions.push_back(torsion);
                }
                for (bond torsion: this->terminal_torsions){
                    this->torsions.push_back(torsion);
                }
                for (bond torsion: this->methylalike_torsions){
                    this->torsions.push_back(torsion);
                }
                valid_mode_input = true;
                break;
            case 2:
                for (bond torsion: this->central_torsions){
                    this->torsions.push_back(torsion);
                }
                valid_mode_input = true;
                break;
            case 3:
                for (bond torsion: this->central_torsions){
                    this->torsions.push_back(torsion);
                }
                for (bond torsion: this->terminal_torsions){
                    this->torsions.push_back(torsion);
                }
                valid_mode_input = true;
                break;
            default:
                std::cin.clear();
                std::cin.ignore();
                std::cout << "Invalid input." << std::endl;
                break;
        }
    }
    bool valid_increment = false;
    int increment;
    bool valid_start_input = false;
    int input_start;
    bool valid_confirmation = false;
    int confirmation;
    int n_possible_conformers;
    while (!valid_start_input){
        while (!valid_increment){
            std::cout << "Type in an angle increment (30, 45, 60, 90, 120 or 180 in degrees): ";
            std::cin >> increment;
            switch (increment){
                case 30:
                    this->angle_increments.push_back(30.0);
                    this->angle_increments.push_back(45.0);
                    valid_increment = true;
                    break;
                case 45:
                    this->angle_increments.push_back(45.0);
                    this->angle_increments.push_back(60.0);
                    valid_increment = true;
                    break;
                case 60:
                    this->angle_increments.push_back(60.0);
                    this->angle_increments.push_back(90.0);
                    valid_increment = true;
                    break;
                case 90:
                    this->angle_increments.push_back(90.0);
                    this->angle_increments.push_back(120.0);
                    valid_increment = true;
                    break;
                case 120:
                    this->angle_increments.push_back(120.0);
                    this->angle_increments.push_back(180.0);
                    valid_increment = true;
                    break;
                case 180:
                    this->angle_increments.push_back(180.0);
                    valid_increment = true;
                    break;
                default:
                    std::cin.clear();
                    std::cin.ignore();
                    std::cout << "Invalid input." << std::endl;
                    break;
            }
        }
        n_possible_conformers = 0;
        for (double increment: this->angle_increments){
            n_possible_conformers = n_possible_conformers + pow(360/(int)increment, this->torsions.size());

        }
        valid_confirmation = false;
        while (!valid_confirmation){
            std::cout << "Up to " << n_possible_conformers << " will be generated. Start calculation (1) or new increment selection (2)?: ";
            std::cin >> confirmation;
            switch (confirmation){
                case 1:
                    valid_confirmation = true;
                    valid_start_input = true;
                    break;
                case 2:
                    valid_confirmation = true;
                    valid_increment = false;
                    break;
                default:
                    std::cin.clear();
                    std::cin.ignore();
                    std::cout << "Invalid input." << std::endl;
                    break;
            }
        }
    }
    // *** END OF FUNCTION SELECTION MENU ***

    this->generation_setup();

    int i;
    int n_generated_conformers = 0;
    for (double increment: this->angle_increments){
        this->angles.clear();
        for (i = 0; i < 360/(int)increment; i++){
            this->angles.push_back((double)i*increment);
        }
        n_generated_conformers = this->combinations(this->input_coords, 0, n_generated_conformers);
    }

    std::cout << n_generated_conformers << " have been generated." << std::endl;

    return;
}


void ConformerGenerator::get_torsions(){
    int i, j;
    std::string element_i, element_j;

    for (bond bond: this->mol->bonds){
        if (bond.bond_order != 1){
            continue;
        }
        i = bond.atom_index1;
        j = bond.atom_index2;
        element_i = this->mol->atoms[i]->element;
        element_j = this->mol->atoms[j]->element;
        if (is_terminal_atom(element_i) || is_terminal_atom(element_j)){
            continue;
        }
        if (this->mol->atoms[i]->core_of_terminal_group){
            if (element_i == "C"){
                this->methylalike_torsions.push_back(bond);
            }
            else{
                this->terminal_torsions.push_back(bond);
            }
        }
        else if (this->mol->atoms[j]->core_of_terminal_group){
            if (element_j == "C"){
                this->methylalike_torsions.push_back(bond);
            }
            else{
                this->terminal_torsions.push_back(bond);
            }
        }
        else{
            this->central_torsions.push_back(bond);
        }
    }
    
    return;
}


void ConformerGenerator::find_cycles(){
    int ancestors[this->mol->atoms.size()];
    memset(ancestors, -1, sizeof(ancestors));
    char status[this->mol->atoms.size()];
    memset(status, 'U', sizeof(status));
    this->cycle_detection(0, 0, status, ancestors);
    return;
}


// U = UNKNOWN
// V = VISITED
// K = KNOWN
void ConformerGenerator::cycle_detection(int current, int last, char status[], int ancestors[]){
    if (status[current] == 'K' || is_terminal_atom(this->mol->atoms[current]->element)){
        return;
    }
    else if (status[current] == 'V'){
        int i;
        int temp, before;
        int atom1, atom2;
        for (i = 0; i < this->central_torsions.size(); i++){
            atom1 = this->central_torsions[i].atom_index1;
            atom2 = this->central_torsions[i].atom_index2;
            if ((current == atom1 && last == atom2) || (last == atom1 && current == atom2)){
                this->central_torsions.erase(this->central_torsions.begin() + i);
                break;
            }
        }
        temp = last;
        while (temp != current){
            before = temp;
            temp = ancestors[temp];
            for (i = 0; i < this->central_torsions.size(); i++){
                atom1 = this->central_torsions[i].atom_index1;
                atom2 = this->central_torsions[i].atom_index2;
                if ((temp == atom1 && before == atom2) || (before == atom1 && temp == atom2)){
                    this->central_torsions.erase(this->central_torsions.begin() + i);
                    break;
            }
        }
        }
    }
    else if (status[current] == 'U'){
        ancestors[current] = last;
        status[current] = 'V';
        for (int bond_partner: this->mol->atoms[current]->bond_partners){
            if (!is_terminal_atom(this->mol->atoms[bond_partner]->element)){
                if (bond_partner != ancestors[current]){
                    this->cycle_detection(bond_partner, current, status, ancestors);
                }
            }
        }
        status[current] = 'K';
        return;
    }
}


void ConformerGenerator::find_peptidebonds(){
    int i, j, n;
    int peptide1, peptide2;
    int atom1, atom2, atom3;
    int input;
    std::vector<double> coords1;
    std::vector<double> coords2;
    std::vector<double> coords3;
    float distance;
    float distance_double_bond = valence_radii_double["C"] + valence_radii_double["O"] + 0.08;
    bool valid_input;
    std::string element1, element2, element3;
    std::vector<std::vector<int>> peptitebonds;

    for (i = 0; i < this->central_torsions.size(); i++){
        atom1 = this->central_torsions[i].atom_index1;
        atom2 = this->central_torsions[i].atom_index2;
        element1 = this->mol->atoms[atom1]->element;
        element2 = this->mol->atoms[atom2]->element;
        if (element1 == "C" && element2 == "N"){
            for (int atom3: this->mol->atoms[atom1]->bond_partners){
                element3 = this->mol->atoms[atom3]->element;
                if (element3 == "O"){
                    coords1 = this->mol->atoms[atom1]->coords;
                    coords3 = this->mol->atoms[atom3]->coords;
                    distance = sqrt(
                                pow(coords1[0]- coords3[0], 2) +
                                pow(coords1[1]- coords3[1], 2) +
                                pow(coords1[2]- coords3[2], 2)
                               );
                    if (distance <= distance_double_bond){
                        peptitebonds.push_back({atom1, atom3});
                    }
                }
            }
        }
        else if (element1 == "N" && element2 == "C"){
            for (int atom3: this->mol->atoms[atom2]->bond_partners){
                element3 = this->mol->atoms[atom3]->element;
                if (element3 == "O"){
                    coords2 = this->mol->atoms[atom2]->coords;
                    coords3 = this->mol->atoms[atom3]->coords;
                    distance = sqrt(
                                pow(coords2[0]- coords3[0], 2) +
                                pow(coords2[1]- coords3[1], 2) +
                                pow(coords2[2]- coords3[2], 2)
                               );
                    if (distance <= distance_double_bond){
                        peptitebonds.push_back({atom2, atom3});
                    }
                }
            }
        }
    }

    if (!peptitebonds.empty()){
        valid_input = false;
        while (!valid_input){
            std::cout << "Consider peptidebonds (1) or ignore them (2)?: ";
            std::cin >> input;
            if (input == 1){
                valid_input = true;
            }
            else if (input == 2){
                for (i = 0; i < peptitebonds.size(); i++){
                    peptide1 = peptitebonds[i][0];
                    peptide2 = peptitebonds[i][1];
                    n = this->central_torsions.size();
                    for (j = 0; j < n; j++){
                        atom1 = this->central_torsions[j].atom_index1;
                        atom2 = this->central_torsions[j].atom_index2;
                        if ((peptide1 == atom1 && peptide2 == atom2) || (peptide1 == atom2 && peptide2 == atom1)){
                            this->central_torsions.erase(this->central_torsions.begin() + j);
                            break;
                        }
                    }
                }
                valid_input = true;
            }
            else{
                std::cin.clear();
                std::cin.ignore();
                std::cout << "Invalid input." << std::endl;
            }
        }
    }

    return;
}


void ConformerGenerator::generation_setup(){
    int i, j;
    int atom1, atom2;
    int status[this->mol->n_atoms];
    std::vector<int> left_atoms;
    std::vector<int> right_atoms;

    for (bond torsion: this->torsions){
        atom1 = torsion.atom_index1;
        atom2 = torsion.atom_index2;
        memset(status, 0, sizeof(status));
        left_atoms.clear();
        left_atoms = this->torsion_atom_counter(atom1, atom2, status, left_atoms);
        memset(status, 0, sizeof(status));
        right_atoms.clear();
        right_atoms = this->torsion_atom_counter(atom2, atom1, status, right_atoms);
        if (left_atoms.size() <= right_atoms.size()){
            this->torsion_atoms.push_back(std::vector<int>());
            i = this->torsion_atoms.size() - 1;
            for (int j: left_atoms){
                this->torsion_atoms[i].push_back(j);
            }
        }
        else{
            this->torsion_atoms.push_back(std::vector<int>());
            i = this->torsion_atoms.size() - 1;
            for (int j: right_atoms){
                this->torsion_atoms[i].push_back(j);
            }
        }
    }

    for (std::shared_ptr<atom> atom_ptr: this->mol->atoms){
        Eigen::Vector3d coords(atom_ptr->coords.data());
        this->input_coords.push_back(coords);
    }

    return;
}


std::vector<int> ConformerGenerator::torsion_atom_counter(int start, int last, int* status, std::vector<int> container){
    if (status[start] == 1){
        return container;
    }
    else if (is_terminal_atom(this->mol->atoms[start]->element)){
        status[start] = 1;
        container.push_back(start);
        return container;
    }
    else{
        status[start] = 1;
        container.push_back(start);
        for (int bond_partner: this->mol->atoms[start]->bond_partners){
            if (bond_partner != last){
                container = this->torsion_atom_counter(bond_partner, start, status, container);
            }
        }
    }
    return container;
}


int ConformerGenerator::combinations(std::vector<Eigen::Vector3d> new_coords, int index, int counter){
    if (index == this->torsions.size()){
        if (!this->clashes(new_coords)){
            std::string current_workdir = this->workdir_name + std::to_string(counter) + "/";
            std::string coord_file = current_workdir + "coord";
            std::string control_file = current_workdir + "control";
            std::string new_struc = current_workdir + this->struc_filename;
            // create new working directory
            system(("mkdir " + current_workdir).c_str());
            // write new conformer structure to be optimized
            std::ofstream file;
            file.open(new_struc);
            file << this->mol->n_atoms;
            file << "\n\n";
            for (int i = 0; i < this->mol->n_atoms; i++){
                file << this->mol->atoms[i]->element << "   " 
                     << new_coords[i][0]             << "   " 
                     << new_coords[i][1]             << "   " 
                     << new_coords[i][2]             << "\n";
            }
            file.close();
            // convert .xyz file of structure to be optimized to coord file
            system(("x2t " + new_struc + " > " + coord_file).c_str());
            // write control file for UFF optimization
            file.open(control_file);
            file << "$symmetry c1\n";
            file << "$uff\n";
            file << "      2500         1          0 ! maxcycle,modus,nqeq\n";
            file << "    111111                      ! iterm\n";
            file << "  0.10D-07  0.10D-04            ! econv,gconv\n";
            file << "      0.00  1.10                ! qtot,dfac\n";
            file << "  0.10D+03  0.10D-04       0.30 ! epssteep,epssearch,dqmax\n";
            file << "        25      0.10       0.00 ! mxls,dhls,ahls\n";
            file << "      1.00      0.00       0.00 ! alpha,beta,gamma\n";
            file << "         F         F          F ! transform,lnumhess,lmd\n";
            file << "$end\n";
            file.close();
            // perform UFF optimization
            system(("cd " + current_workdir + " ; uff > uff.out 2>&1 &").c_str());
            // cout new conformer
            return counter+1;
        }
        else{
            return counter;
        }
    }
    else{
        int atom1 = this->torsions[index].atom_index1;
        int atom2 = this->torsions[index].atom_index2;
        Eigen::Vector3d axis_vec1(new_coords[atom1].data());
        Eigen::Vector3d axis_vec2(new_coords[atom2].data());
        Eigen::Vector3d axis = axis_vec2 - axis_vec1;
        axis.normalize();
        Eigen::Vector3d new_coord;
        for (double deg: this->angles){
            double rad = 2*M_PI*deg/360.0;
            std::vector<Eigen::Vector3d> new_coords_copy = new_coords;
            for (int torsion_atom: this->torsion_atoms[index]){
                new_coord = new_coords_copy[torsion_atom];
                new_coord = new_coord - axis_vec1;
                new_coord = axis.dot(new_coord) * axis
                            + cos(rad) * axis.cross(new_coord).cross(axis) // FUNKTIONIERT cross().cross()??
                            + sin(rad) * axis.cross(new_coord);
                new_coord = new_coord + axis_vec1;
                new_coords_copy[torsion_atom] = new_coord;
            }
            counter = this->combinations(new_coords_copy, index+1, counter);
        }
        return counter;
    }
}


void ConformerGenerator::write_xyz(std::vector<Eigen::Vector3d> coords, std::string folder, int structure_number){
    int i = 0;
    std::string file = "conformer" + std::to_string(structure_number) + ".xyz";
    file = folder + file;
    std::ofstream new_xyz_file(file);
    new_xyz_file << this->mol->n_atoms;
    new_xyz_file << "\n\n";
    for (Eigen::Vector3d coord: coords){
        new_xyz_file << this->mol->atoms[i]->element << "   " 
                     << coord[0] << "   " 
                     << coord[1] << "   " 
                     << coord[2] << "\n";
        i++;
    }
    new_xyz_file.close();
}


bool ConformerGenerator::clashes(std::vector<Eigen::Vector3d> coords){
    int i, j;
    std::string element_i, element_j;
    double dist, min_dist;
    double tolerance = 0.08;
    
    for (i = 0; i < coords.size()-1; i++){
        for (j = i + 1; j < coords.size(); j++){
            element_i = this->mol->atoms[i]->element;
            element_j = this->mol->atoms[j]->element;
            dist = (coords[i] - coords[j]).norm();
            min_dist = valence_radii_single[element_i] + valence_radii_single[element_j] + tolerance;
            if (dist < min_dist){
                if (this->distant_atoms(i, j)){
                    return true;
                }
            }
        }
    }

    return false;
}


bool ConformerGenerator::distant_atoms(int atom1, int atom2){
    int dist = 1;
    int n_neighbors;
    int curr;
    std::queue<int> neighbors;

    for (int bond_partner: this->mol->atoms[atom1]->bond_partners){
        neighbors.push(bond_partner);
    }
    while (dist < 4 && !neighbors.empty()){
        n_neighbors = neighbors.size();
        while (n_neighbors){
            curr = neighbors.front();
            neighbors.pop();
            if (curr == atom2){
                return false;
            }
            else{
                for (int bond_partner: this->mol->atoms[curr]->bond_partners){
                    neighbors.push(bond_partner);
                }
            }
            n_neighbors--;
        }
        dist++;
    }
    
    return true;
}