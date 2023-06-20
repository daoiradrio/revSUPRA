#include <conformergenerator.hpp>



ConformerGenerator::ConformerGenerator(Structure input_mol){
    char buffer[BUFFER_SIZE];
    this->curr_work_dir = getcwd(buffer, BUFFER_SIZE);
    this->curr_work_dir = this->curr_work_dir + "/";
    this->mol = std::make_shared<Structure>(input_mol);
};



ConformerGenerator::~ConformerGenerator(){};



void ConformerGenerator::generate_conformers(){
    this->get_torsions();

    this->find_cycles();

    this->find_peptidebonds();

    this->selection_menu();

    this->generation_setup();

    std::cout << "Generating conformer structures..." << std::endl;
    int i;
    int n_generated_conformers = 0;
    for (int increment: this->angle_increments){
        this->angles.clear();
        for (i = 0; i < 360/increment; i++){
            this->angles.push_back(i*increment);
        }
        this->check_rot_sym(increment);
	    //n_generated_conformers = this->combinations(this->input_coords, 0, n_generated_conformers);
	    //n_generated_conformers = this->old_combinations(this->input_coords_mat, 0, n_generated_conformers);
        n_generated_conformers = this->combinations(this->input_coords_mat, 0, n_generated_conformers);
    }
    std::string command;
    if (n_generated_conformers){
        command = "mkdir " + this->output_foldername;
        system(command.c_str());
        command = "mv " + this->curr_work_dir + this->struc_filename + "* " + this->curr_work_dir + this->output_foldername;
        system(command.c_str());
    }
    this->analyzer.remove_doubles(this->output_foldername, this->struc_filename);

    return;
}



void ConformerGenerator::get_torsions()
{
    std::string                 element_1, element_2;
    std::shared_ptr<Torsion>    new_torsion;

    for (const std::shared_ptr<Bond>& bond: this->mol->bonds){
        if (bond->bond_order != 1){
            continue;
        }
        element_1 = bond->atom1->element;
        element_2 = bond->atom2->element;
        if (is_terminal_atom(element_1) || is_terminal_atom(element_2)){
            continue;
        }
        new_torsion = std::make_shared<Torsion>();
        //new_torsion->bond = std::move(bond); // AUSPROBIEREN
        new_torsion->bond = bond;
        if (bond->atom1->core_of_terminal_group){
            if (element_1 == "C"){
                this->methylalike_torsions.push_back(std::move(new_torsion));
            }
            else{
                this->terminal_torsions.push_back(std::move(new_torsion));
            }
        }
        else if (bond->atom2->core_of_terminal_group){
            if (element_2 == "C"){
                this->methylalike_torsions.push_back(std::move(new_torsion));
            }
            else{
                this->terminal_torsions.push_back(std::move(new_torsion));
            }
        }
        else{
            this->central_torsions.push_back(std::move(new_torsion));
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
            atom1 = this->central_torsions[i]->bond->atom1->index;
            atom2 = this->central_torsions[i]->bond->atom2->index;
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
                atom1 = this->central_torsions[i]->bond->atom1->index;
                atom2 = this->central_torsions[i]->bond->atom2->index;
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
        atom1 = this->central_torsions[i]->bond->atom1->index;
        atom2 = this->central_torsions[i]->bond->atom2->index;
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
                        atom1 = this->central_torsions[j]->bond->atom1->index;
                        atom2 = this->central_torsions[j]->bond->atom2->index;
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



void ConformerGenerator::selection_menu(){
    bool start_calculation;
    bool confirm_increment;
    bool valid_mode_input;
    bool valid_increment;
    bool flag_n_conformers;

    int mode_input;
    int increment_input;
    int confirm_input;

    int n_possible_conformers;

    start_calculation = false;
    while (!start_calculation){
        valid_mode_input = false;
	this->torsions.clear();
        while (!valid_mode_input){
            std::cout << "Consider rotatable all bonds to terminal groups like -CH3, -NH2, -OH (1) "
                         "or ignore them (2) "
                         "or ignore just ignore bonds to methyl-like groups (3)?: ";
            std::cin >> mode_input;
            switch (mode_input){
                case 1:
                    for (std::shared_ptr<Torsion>& torsion: this->central_torsions){
                        this->torsions.push_back(torsion);
                    }
                    for (std::shared_ptr<Torsion>& torsion: this->terminal_torsions){
                        this->torsions.push_back(torsion);
                    }
                    for (std::shared_ptr<Torsion>& torsion: this->methylalike_torsions){
                        this->torsions.push_back(torsion);
                    }
                    valid_mode_input = true;
                    break;
                case 2:
                    for (std::shared_ptr<Torsion>& torsion: this->central_torsions){
                        this->torsions.push_back(torsion);
                    }
                    valid_mode_input = true;
                    break;
                case 3:
                    for (std::shared_ptr<Torsion>& torsion: this->central_torsions){
                        this->torsions.push_back(torsion);
                    }
                    for (std::shared_ptr<Torsion>& torsion: this->terminal_torsions){
                        this->torsions.push_back(torsion);
                    }
                    valid_mode_input = true;
                    break;
                default:
                    std::cin.clear();
                    std::cin.ignore();
                    std::cout << "Invalid input." << std::endl;
            }
        }
        confirm_increment = false;
        this->angle_increments.clear();
	    valid_increment = false;
        while (!confirm_increment){
            std::cout << "Type in an angle increment (30, 45, 60, 90, 120 or 180 in degrees): ";
            std::cin >> increment_input;
            switch (increment_input){
                case 30:
                    this->angle_increments.push_back(30);
                    this->angle_increments.push_back(45);
                    valid_increment = true;
                    break;
                case 45:
                    this->angle_increments.push_back(45);
                    this->angle_increments.push_back(60);
                    valid_increment = true;
                    break;
                case 60:
                    this->angle_increments.push_back(60);
                    this->angle_increments.push_back(90);
                    valid_increment = true;
                    break;
                case 90:
                    this->angle_increments.push_back(90);
                    this->angle_increments.push_back(120);
                    valid_increment = true;
                    break;
                case 120:
                    this->angle_increments.push_back(120);
                    this->angle_increments.push_back(180);
                    valid_increment = true;
                    break;
                case 180:
                    this->angle_increments.push_back(180);
                    valid_increment = true;
                    break;
                default:
                    std::cin.clear();
                    std::cin.ignore();
                    std::cout << "Invalid input." << std::endl;
            }
            if (valid_increment){
                n_possible_conformers = 0;
                for (int increment: this->angle_increments){
                    n_possible_conformers = n_possible_conformers + pow(360/increment, this->torsions.size());
                }
                flag_n_conformers = false;
                while (!flag_n_conformers){
                    std::cout << "Up to " << n_possible_conformers << " will be generated. "
                                 "Start calculation (1) or restart selection (2)? ";
                    std::cin >> confirm_input;
                    switch (confirm_input){
                        case 1:
                            flag_n_conformers = true;
                            confirm_increment = true;
                            start_calculation = true;
                            break;
                        case 2:
                            flag_n_conformers = true;
                            confirm_increment = true;
                        default:
                            std::cin.clear();
                            std::cin.ignore();
                            std::cout << "Invalid input." << std::endl;
                            break;
                    }
                }
            }
        }
    }
}



void ConformerGenerator::generation_setup(){
    int 		        i, j;
    int 		        atom1, atom2;
    std::vector<int>	left_atoms;
    std::vector<int>	right_atoms;
    std::vector<int>	status(this->mol->n_atoms, 0);
    Symmetry 		    sym;
    std::vector<int> 	torsion_group_left;
    std::vector<int> 	torsion_group_right;

    this->optimizer = Optimizer();
    this->analyzer  = Analyzer();

    this->input_coords_mat = this->mol->coords;

    for (i = 0; i < this->mol->coords.rows(); i++){
	    this->input_coords.push_back(this->mol->coords.row(i));
    }

    for (std::shared_ptr<Torsion>& torsion: this->torsions){
        atom1 = torsion->bond->atom1->index;
        atom2 = torsion->bond->atom2->index;
        left_atoms.clear();
        std::fill(status.begin(), status.end(), 0);
        left_atoms = this->torsion_atom_counter(atom1, atom2, status, left_atoms);
        torsion->rot_atoms1 = left_atoms;
        right_atoms.clear();
        std::fill(status.begin(), status.end(), 0);
        right_atoms = this->torsion_atom_counter(atom2, atom1, status, right_atoms);
        torsion->rot_atoms2 = right_atoms;
        if (left_atoms.size() <= right_atoms.size()){
            torsion->rot_atoms = left_atoms;
        }
        else{
            torsion->rot_atoms = right_atoms;
        }
    }

    return;
}



std::vector<int> ConformerGenerator::torsion_atom_counter(int start, int last, std::vector<int> status, std::vector<int> container){
    if (status[start] == 1){
        std::sort(container.begin(), container.end());
        return container;
    }
    else{
        status[start] = 1;
        for (int bond_partner: this->mol->atoms[start]->bond_partners){
            if (bond_partner != last){
                container.push_back(bond_partner);
                container = this->torsion_atom_counter(bond_partner, start, status, container);
            }
        }
    }
    std::sort(container.begin(), container.end());
    return container;
}



std::vector<int> ConformerGenerator::get_torsion_group(int start, int last, std::vector<int> status, std::vector<int> container){
    if (status[start] == 1){
        std::sort(container.begin(), container.end());
        return container;
    }
    else{
        status[start] = 1;
	    for (int bond_partner: this->mol->atoms[start]->bond_partners){
	        if (bond_partner != last){
		        if (std::find_if(
		            this->torsions.begin(), this->torsions.end(),
			        [&](const std::shared_ptr<Torsion>& torsion){
			            return (torsion->bond->atom1->index == bond_partner || torsion->bond->atom2->index == bond_partner);
			        }
		        ) == this->torsions.end()){
	                container.push_back(bond_partner);
		            container = this->get_torsion_group(bond_partner, start, status, container);
	            }
	        }
	    }
    }
    std::sort(container.begin(), container.end());
    return container;
}



void ConformerGenerator::check_rot_sym(int angle_increment)
{
    Symmetry                    sym;
    int                         i, j, k;
    std::vector<int>            left_torsion_group;
    std::vector<int>            right_torsion_group;
    std::vector<int>	        status(this->mol->n_atoms, 0);
    std::vector<int>            torsion_done(this->torsions.size(), 0);
    std::shared_ptr<Torsion>    torsion1;
    std::shared_ptr<Torsion>    torsion2;

    for (std::shared_ptr<Torsion>& torsion: this->torsions){
        std::fill(status.begin(), status.end(), 0);
        left_torsion_group = this->get_torsion_group(
            torsion->bond->atom1->index,
            torsion->bond->atom2->index,
            status
        );
        torsion->rot_sym1 = sym.rot_sym_along_bond(
            this->mol,
            left_torsion_group,
            torsion->bond->atom1->index,
            torsion->bond->atom2->index
        );
        if (torsion->rot_sym1 > 1){
            torsion->rot_sym_atoms1 = left_torsion_group;
        }
        std::fill(status.begin(), status.end(), 0);
        right_torsion_group = this->get_torsion_group(
            torsion->bond->atom2->index,
            torsion->bond->atom1->index,
            status
        );
        torsion->rot_sym2 = sym.rot_sym_along_bond(
            this->mol,
            right_torsion_group,
            torsion->bond->atom2->index,
            torsion->bond->atom1->index
        );
        if (torsion->rot_sym2 > 1){
            torsion->rot_sym_atoms2 = right_torsion_group;
        }
    }

    this->rot_angles.clear();
    for (const std::shared_ptr<Torsion>& torsion: this->torsions){
        this->rot_angles.push_back({});
    }

    std::fill(torsion_done.begin(), torsion_done.end(), 0);
    
    for (i = 0; i < this->torsions.size(); i++){
        if (torsion_done[i]){
            continue;
        }
        torsion1 = this->torsions[i];
        if (torsion1->rot_sym1 == 1 || std::max(360/torsion1->rot_sym1, angle_increment) % std::min(360/torsion1->rot_sym1, angle_increment) != 0){
            if (torsion1->rot_sym2 == 1 || std::max(360/torsion1->rot_sym2, angle_increment) % std::min(360/torsion1->rot_sym2, angle_increment) != 0){
                for (k = 0; k < 360/angle_increment; k++){
                    this->rot_angles[i].push_back(k*angle_increment);
                }
                torsion_done[i] = 1;
                continue;
            }
        }
        if (torsion1->rot_sym_atoms1 == torsion1->rot_atoms1){
            if (torsion1->rot_sym1 > 1){
                if (std::max(360/torsion1->rot_sym1, angle_increment)%std::min(360/torsion1->rot_sym1, angle_increment) == 0){
                    for (j = 0; j*angle_increment < 360/torsion1->rot_sym1; j++){
                        this->rot_angles[i].push_back(j*angle_increment);
                    }
                    torsion_done[i] = 1;
                    continue;
                }
            }
        }
        if (torsion1->rot_sym_atoms2 == torsion1->rot_atoms2){
            if (torsion1->rot_sym2 > 1){
                if (std::max(360/torsion1->rot_sym2, angle_increment)%std::min(360/torsion1->rot_sym2, angle_increment) == 0){
                    for (j = 0; j*angle_increment < 360/torsion1->rot_sym2; j++){
                        this->rot_angles[i].push_back(j*angle_increment);
                    }
                    torsion_done[i] = 1;
                    continue;
                }
            }
        }
        for (j = i+1; j < this->torsions.size(); j++){
            torsion2 = this->torsions[j];
            if (torsion1->rot_sym1 > 1){
                if (std::max(360/torsion1->rot_sym1, angle_increment)%std::min(360/torsion1->rot_sym1, angle_increment) == 0){
                    if (torsion1->rot_sym_atoms1 == torsion2->rot_sym_atoms1 || torsion1->rot_sym_atoms1 == torsion2->rot_sym_atoms2){
                        if (!torsion_done[j]){
                            for (k = 0; k*angle_increment < 360/torsion1->rot_sym1; k++){
                                this->rot_angles[j].push_back(k*angle_increment);
                            }
                            torsion_done[j] = 1;
                            for (k = 0; k < 360/angle_increment; k++){
                                this->rot_angles[i].push_back(k*angle_increment);
                            }
                            torsion_done[i] = 1;
                            continue;
                        }
                    }
                }
            }
            if (torsion1->rot_sym2 > 1){
                if (std::max(360/torsion1->rot_sym2, angle_increment)%std::min(360/torsion1->rot_sym2, angle_increment) == 0){
                    if (torsion1->rot_sym_atoms2 == torsion2->rot_sym_atoms1 || torsion1->rot_sym_atoms2 == torsion2->rot_sym_atoms2){
                        if (!torsion_done[j]){
                            for (k = 0; k*angle_increment < 360/torsion1->rot_sym2; k++){
                                this->rot_angles[j].push_back(k*angle_increment);
                            }
                            torsion_done[j] = 1;
                            for (k = 0; k < 360/angle_increment; k++){
                                this->rot_angles[i].push_back(k*angle_increment);
                            }
                            torsion_done[i] = 1;
                            continue;
                        }
                    }
                }
            }
        }
    }
    std::cout << "Inkrement " << angle_increment << std::endl;
    std::cout << std::endl;
    i = 0;
    for (std::shared_ptr<Torsion> torsion: this->torsions){
         std::cout << torsion->bond->atom1->element << torsion->bond->atom1->index << " "
                   << torsion->bond->atom2->element << torsion->bond->atom2->index << std::endl;
        for (int angle: this->rot_angles[i]){
            std::cout << angle << " ";
        }
        std::cout << "\n\n";
        i++;
    }

    return;
}



int ConformerGenerator::old_combinations(Eigen::MatrixX3d new_coords, int index, int counter){
    if (index == this->torsions.size()){
        if (!this->clashes(new_coords)){
            // get necessary file paths for working directory of new conformer structure
	        int fin;
            std::string current_workdir = this->workdir_name + std::to_string(counter) + "/";
            std::string coord_file = current_workdir + "coord";
            std::string control_file = current_workdir + "control";
            std::string new_struc = current_workdir + this->struc_filename + std::to_string(counter);
	        std::string opt_struc = current_workdir + this->opt_struc_filename;
	        std::string command;
            // create new working directory
	        command = "mkdir " + current_workdir;
	        system(command.c_str());
            // write new conformer structure to be optimized
            this->write_xyz(new_coords, new_struc);
            // convert .xyz file of structure to be optimized to coord file
	        command = "x2t " + new_struc + " > " + coord_file;
            system(command.c_str());
            // write control file for UFF optimization
	        std::ofstream file;
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
	        command = "cd " + current_workdir + " ; uff > uff.out 2>&1 &";
            fin = system(command.c_str());
            command = "t2x " + coord_file + " > " + new_struc + " 2>/dev/null";
            fin = system(command.c_str());
            command = "mv ~/revSUPRA/" + new_struc + " ~/revSUPRA/conformer" + std::to_string(counter) + ".xyz";
            fin = system(command.c_str());
            command = "rm -rf " + current_workdir;
            fin = system(command.c_str());
            return counter+1;
        }
        else{
            return counter;
        }
        return counter;
    }
    else{
        int atom1 = this->torsions[index]->bond->atom1->index;
        int atom2 = this->torsions[index]->bond->atom2->index;
        for (int angle: this->angles){
	        RotationAxis rot_axis(new_coords.row(atom1), new_coords.row(atom2));
            Eigen::MatrixX3d new_coords_copy = new_coords;
            Eigen::Vector3d new_coord;
            for (int torsion_atom: this->torsions[index]->rot_atoms){
            //for (int torsion_atom: this->torsion_atoms[index]){
                new_coord = rot_axis.rotate_atom(new_coords_copy.row(torsion_atom), angle);
                new_coords_copy.row(torsion_atom) = new_coord;
            }
            counter = this->old_combinations(new_coords_copy, index+1, counter);
        }
        return counter;
    }
}



int ConformerGenerator::combinations(std::vector<Eigen::Vector3d> new_coords, int index, int counter){
    if (index == this->torsions.size()){
        if (!this->clashes(new_coords)){
            // get necessary file paths for working directory of new conformer structure
	        int fin;
            std::string current_workdir = this->workdir_name + std::to_string(counter) + "/";
            std::string coord_file = current_workdir + "coord";
            std::string control_file = current_workdir + "control";
            std::string new_struc = current_workdir + this->struc_filename;
	        std::string opt_struc = current_workdir + this->opt_struc_filename;
	        std::string command;
            // create new working directory
	        command = "mkdir " + current_workdir;
	        system(command.c_str());
            // write new conformer structure to be optimized
            this->write_xyz(new_coords, new_struc);
            // convert .xyz file of structure to be optimized to coord file
	        command = "x2t " + new_struc + " > " + coord_file;
            system(command.c_str());
            // write control file for UFF optimization
	        std::ofstream file;
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
	        command = "cd " + current_workdir + " ; uff > uff.out 2>&1 &";
            fin = system(command.c_str());
            command = "t2x " + coord_file + " > " + new_struc + " 2>/dev/null";
            fin = system(command.c_str());
            command = "mv ~/revSUPRA/" + new_struc + " ~/revSUPRA/conformer" + std::to_string(counter) + ".xyz";
            fin = system(command.c_str());
            command = "rm -rf " + current_workdir;
            fin = system(command.c_str());
            // cout new conformer
            return counter+1;
        }
        else{
            return counter;
        }
    }
    else{
        int atom1 = this->torsions[index]->bond->atom1->index;
        int atom2 = this->torsions[index]->bond->atom2->index;
        //Eigen::Vector3d axis_vec1(new_coords[atom1].data());
        //Eigen::Vector3d axis_vec2(new_coords[atom2].data());
	    Eigen::Vector3d axis_vec1(new_coords[atom1]);
	    Eigen::Vector3d axis_vec2(new_coords[atom2]);
        Eigen::Vector3d axis = axis_vec2 - axis_vec1;
        axis.normalize();
        Eigen::Vector3d new_coord;
        for (int deg: this->angles){
            double rad = 2*M_PI*(double)deg/360.0;
            std::vector<Eigen::Vector3d> new_coords_copy = new_coords;
            for (int torsion_atom: this->torsions[index]->rot_atoms){
            //for (int torsion_atom: this->torsion_atoms[index]){
                new_coord = new_coords_copy[torsion_atom];
                new_coord = new_coord - axis_vec1;
                new_coord = axis.dot(new_coord) * axis
                            + cos(rad) * axis.cross(new_coord).cross(axis)
                            + sin(rad) * axis.cross(new_coord);
                new_coord = new_coord + axis_vec1;
                new_coords_copy[torsion_atom] = new_coord;
            }
            counter = this->combinations(new_coords_copy, index+1, counter);
        }
        return counter;
    }
}



int ConformerGenerator::combinations(Eigen::MatrixX3d new_coords, int index, int counter){
    if (index == this->torsions.size()){
        if (!this->clashes(new_coords)){
            std::string new_struc = this->struc_filename + std::to_string(counter) + ".xyz";
            this->write_xyz(new_coords, new_struc);
            int fin = this->optimizer.uff_optimization(this->curr_work_dir, new_struc, counter);
            return counter+1;
        }
        else{
            return counter;
        }
    }
    else{
        int atom1 = this->torsions[index]->bond->atom1->index;
        int atom2 = this->torsions[index]->bond->atom2->index;
        //for (int angle: this->angles){
            //if (index == 3 && (angle == 120 || angle == 240)){
            //    continue;
            //}
        for (int angle: this->rot_angles[index]){
            Eigen::MatrixX3d new_coords_copy = new_coords;
            Eigen::Vector3d new_coord;
            for (int torsion_atom: this->torsions[index]->rot_atoms){
                new_coord = RotationAxis::rotate_atom(
                    new_coords.row(atom1),
                    new_coords.row(atom2),
                    new_coords_copy.row(torsion_atom),
                    angle
                );
                new_coords_copy.row(torsion_atom) = new_coord;
            }
            counter = this->combinations(new_coords_copy, index+1, counter);
        }
        return counter;
    }
}



void ConformerGenerator::write_xyz(std::vector<Eigen::Vector3d> coords, std::string destination){
    std::ofstream file(destination);
    file << this->mol->n_atoms;
    file << "\n\n";
    for (int i = 0; i < this->mol->n_atoms; i++){
        file << this->mol->atoms[i]->element << "   " 
             << coords[i][0]                 << "   " 
             << coords[i][1]                 << "   " 
             << coords[i][2]                 << "\n";
    }
    file.close();
}



void ConformerGenerator::write_xyz(Eigen::MatrixX3d coords, std::string destination){
    std::ofstream file(destination);
    file << this->mol->n_atoms;
    file << "\n\n";
    for (int i = 0; i < this->mol->n_atoms; i++){
        file << this->mol->atoms[i]->element << "   " 
             << coords(i, 0)                 << "   " 
             << coords(i, 1)                 << "   " 
             << coords(i, 2)                 << "\n";
    }
    file.close();
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



bool ConformerGenerator::clashes(Eigen::MatrixX3d coords){
    int i, j;
    std::string element_i, element_j;
    double dist, min_dist;
    double tolerance = 0.08;
    
    for (i = 0; i < coords.rows(); i++){
        for (j = i + 1; j < coords.rows(); j++){
            element_i = this->mol->atoms[i]->element;
            element_j = this->mol->atoms[j]->element;
            dist = (coords.row(i) - coords.row(j)).norm();
            //dist = (coords[i] - coords[j]).norm();
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
