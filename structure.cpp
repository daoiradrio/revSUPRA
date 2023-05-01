#include "structure.hpp"



Structure::Structure(){}


Structure::~Structure(){}


void Structure::get_structure(std::string filepath){
    this->read_xyz(filepath);
    this->get_bond_matrix();
}


void Structure::read_xyz(std::string filepath){
    std::string element;
    std::string new_label;
    std::string line;
    std::ifstream file(filepath);
    atom new_atom;
    double xcoord, ycoord, zcoord;
    int atom_index;
    int line_index;

    if (file.is_open()){
        line_index = 0;
        atom_index = 0;
        while (getline(file, line)){
            std::stringstream linestream(line);
            if (line_index == 0){
                linestream >> this->n_atoms;
            }
            else if (line_index >= 2){
                std::stringstream linestream(line);
                linestream >> element >> xcoord >> ycoord >> zcoord;
                new_atom.element = element;
                new_atom.label = element + std::to_string(atom_index);
                new_atom.coords[0] = xcoord;
                new_atom.coords[1] = ycoord;
                new_atom.coords[2] = zcoord;
                new_atom.bond_partners = {};
                this->atoms.push_back(new_atom);
                atom_index++;
            }
            line_index++;
        }
    }
    else{
        std::cout << "FAILED OPENING .xyz FILE!" << std::endl;
    }
    file.close();

    return;
}


void Structure::get_bond_matrix(){
    int i, j, dim;
    int valence, max_valence;
    int bond_order;
    atom* atom1;
    atom* atom2;

    dim = (this->n_atoms-1)*this->n_atoms/2;
    this->bond_matrix = new int[dim];

    for (i = 0; i < this->n_atoms; i++){
        valence = 0;
        atom1 = &(this->atoms[i]);
        max_valence = max_valences[get_element(atom1->element)];
        for (j = i + 1; j < this->n_atoms; j++){
            if (valence >= max_valence){
                break;
            }
            atom2 = &(this->atoms[j]);
            bond_order = this->get_bond_order(atom1, atom2);
            this->bond_matrix[this->hash_bond_matrix(i, j)] = bond_order;
        }
    }

    return;
}


int Structure::hash_bond_matrix(int row_index, int column_index){
    return (2*row_index*this->n_atoms - row_index - pow(row_index, 2))/2 + column_index - row_index - 1;
}


/*void Structure::get_connectivity(){
    int i, j, n;
    int valence, max_valence;
    int bond_order;
    atom* atom1;
    atom* atom2;
    bond new_bond;
    std::string element1, element2;
    
    for (i = 0; i < this->n_atoms; i++){
        valence = 0;
        atom1 = &(this->atoms[i]);
        max_valence = max_valences[atom1->element];
        for (j = i + 1; j < this->n_atoms; j++){
            if (valence >= max_valence){
                break;
            }
            atom2 = &(this->atoms[j]);
            bond_order = this->get_bond_order(atom1, atom2);
            if (bond_order){
                new_bond.atom1 = *atom1;
                new_bond.atom2 = *atom2;
                new_bond.bond_order = bond_order;
                this->bonds.push_back(new_bond);
                atom1->bond_partners.push_back(*atom2);
                atom2->bond_partners.push_back(*atom1);
            }
        }
    }

    return;
}*/


int Structure::get_bond_order(atom* atom1, atom* atom2)
{
    int bond_order = 0;
    float tolerance = 0.08;
    double single_bond = -1000.0;
    double double_bond = -1000.0;
    double triple_bond = -1000.0;

    float valence_radius_single1 = valence_radii_single[get_element(atom1->label)];
    float valence_radius_single2 = valence_radii_single[get_element(atom2->label)];
    float valence_radius_double1 = valence_radii_double[get_element(atom1->label)];
    float valence_radius_double2 = valence_radii_double[get_element(atom2->label)];
    float valence_radius_triple1 = valence_radii_triple[get_element(atom1->label)];
    float valence_radius_triple2 = valence_radii_triple[get_element(atom2->label)];

    double distance = sqrt(
        pow((atom1->coords[0] - atom2->coords[0]), 2) +
        pow((atom1->coords[1] - atom2->coords[1]), 2) +
        pow((atom1->coords[2] - atom2->coords[2]), 2)
    );

    if (valence_radius_single1 && valence_radius_single2){
        single_bond = valence_radius_single1 + valence_radius_single2 + tolerance;
    }
    if (valence_radius_double1 && valence_radius_double2){
        double_bond = valence_radius_double1 + valence_radius_double2 + tolerance;
    }
    if (valence_radius_triple1 && valence_radius_triple2){
        triple_bond = valence_radius_triple1 + valence_radius_triple2 + tolerance;
    }

    if (distance <= triple_bond){
        bond_order = 3;
    }
    else if (distance <= double_bond){
        bond_order = 2;
    }
    else if (distance <= single_bond){
        bond_order = 1;
    }

    return bond_order;
}