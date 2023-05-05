#include "conformergenerator.hpp"



ConformerGenerator::ConformerGenerator(Structure input_mol){
    this->mol = std::make_shared<Structure>(input_mol);
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
                    break;
            }
        }
        n_possible_conformers = 0;
        for (int increment: this->angle_increments){
            n_possible_conformers = n_possible_conformers + pow(360/increment, this->torsions.size());

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
    for (int increment: this->angle_increments){
        this->angles.clear();
        for (i = 0; i < 360/increment; i++){
            this->angles.push_back(i*increment);
        }

    }

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
    std::vector<float> coords1;
    std::vector<float> coords2;
    std::vector<float> coords3;
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


int ConformerGenerator::combinations(std::vector<std::vector<float>> new_coords, int index, int counter){
    if (index == new_coords.size()){
        return counter+1;
    }
    else{
        int atom1 = this->torsions[index].atom_index1;
        int atom2 = this->torsions[index].atom_index2;
        std::vector<float> axis_vec1 = this->mol->atoms[atom1]->coords;
        std::vector<float> axis_vec2 = this->mol->atoms[atom2]->coords;
        float rad;
        for (int deg: this->angles){
            deg = deg/360;
            rad = 2*M_PI*(float)deg;
        }
        return counter;
    }
}