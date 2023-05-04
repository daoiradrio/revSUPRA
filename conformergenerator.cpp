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