#include "conformergenerator.hpp"



ConformerGenerator::ConformerGenerator(){};


ConformerGenerator::~ConformerGenerator(){};


void ConformerGenerator::generate_conformers(Structure* molecule){
    this->get_torsions(molecule->bonds);
    //this->find_cycles(molecule->atoms);

    return;
}


void ConformerGenerator::get_torsions(std::vector<bond> bonds){
    int i, j, k;
    int bond_order;
    atom atom1;
    atom atom2;
    atom* new1;
    atom* new2;

    /*for (i = 0; i < bonds.size(); i++){
        std::cout << bonds[i].atom1.label << " " << bonds[i].atom2.label << std::endl;
        std::cout << bonds[i].atom1.bond_partners.size() << " " << bonds[i].atom2.bond_partners.size() << std::endl;
    }*/

    /*
    int i, j, k;
    int bond_order;
    atom atom1;
    atom atom2;

    for (bond bond: bonds){
        if (bond.bond_order != 1){
            continue;
        }
        atom1 = bond.atom1;
        atom2 = bond.atom2;
        std::cout << bond.atom1.label << std::endl;
        std::cout << bond.atom1.bond_partners.size() << std::endl;
        std::cout << bond.atom2.label << std::endl;
        std::cout << bond.atom2.bond_partners.size() << std::endl;
        std::cout << "_____________________" << std::endl;
        std::cout << "______________" << std::endl;
        std::cout << atom1->label << std::endl;
        std::cout << atom1->core_of_terminal_group << std::endl;
        std::cout << atom2->label << std::endl;
        std::cout << atom2->core_of_terminal_group << std::endl;
        std::cout << "______________" << std::endl;
        if (is_terminal_atom(atom1->element) || is_terminal_atom(atom2->element)){
            continue;
        }
        if (atom1->core_of_terminal_group){
            if (atom1->element == "C"){
                this->methylalike_torsions.push_back(bond);
            }
            else{
                this->terminal_torsions.push_back(bond);
            }
        }
        else if (atom2->core_of_terminal_group){
            if (atom2->element == "C"){
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
    */
    return;
}

/*
void ConformerGenerator::find_cycles(std::vector<atom*> atoms){
    atom* start = atoms[0];
    int ancestors[atoms.size()];
    memset(ancestors, -1, sizeof(ancestors));
    char status[atoms.size()];
    memset(status, 'U', sizeof(status));
    this->cycle_detection(atoms, start, start, status, ancestors);

    return;
}


// U = UNKNOWN
// V = VISITED
// K = KNOWN
void ConformerGenerator::cycle_detection(
    std::vector<atom*> atoms, atom* current, atom* last, char status[], int ancestors[]
){
    if (status[current->index] == 'K' || is_terminal_atom(current->element)){
        return;
    }
    else if (status[current->index] == 'V'){
        bond* new_bond;
        for (int i = 0; i < this->central_torsions.size(); i++){
            bond* torsion = this->central_torsions[i];
            if (current->label == torsion->atom1->label || current->label == torsion->atom2->label){
                if (last->label == torsion->atom1->label || last->label == torsion->atom2->label){
                    new_bond->atom1 = current;
                    new_bond->atom2 = last;
                    this->ring_bonds.push_back(new_bond);
                    this->central_torsions.erase(this->central_torsions.begin()+i);
                    break;
                }
            }
        }
        atom* temp = last;
        atom* before;
        while (temp->label != current->label){
            before = temp;
            temp = atoms[ancestors[temp->index]];
            for (int i = 0; i < this->central_torsions.size(); i++){
                bond* torsion = this->central_torsions[i];
                if (temp->label == torsion->atom1->label || temp->label == torsion->atom2->label){
                    if (before->label == torsion->atom1->label || before->label == torsion->atom2->label){
                        new_bond->atom1 = temp;
                        new_bond->atom2 = before;
                        this->ring_bonds.push_back(new_bond);
                        this->central_torsions.erase(this->central_torsions.begin()+i);
                        break;
                    }
                }
            }
        }
    }
    else if (status[current->index] == 'U'){
        ancestors[current->index] = last->index;
        status[current->index] = 'V';
        for (atom* bond_partner: current->bond_partners){
            if (!is_terminal_atom(bond_partner->element)){
                if (bond_partner->label != current->label){
                    this->cycle_detection(atoms, bond_partner, current, status, ancestors);
                }
            }
        }
        status[current->index] = 'K';
        return;
    }
}
*/