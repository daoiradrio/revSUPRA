#include "structure.hpp"



Structure::Structure(){}


Structure::~Structure(){}


void Structure::get_structure(std::string filepath){
    this->read_xyz(filepath);
    this->get_bonds();
}


void Structure::read_xyz(std::string filepath){
    std::unique_ptr<atom> new_atom;
    //auto new_atom = std::make_unique<atom>();
    //atom* new_atom;
    std::string element;
    std::string new_label;
    std::string line;
    std::ifstream file(filepath);
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
                new_atom = std::make_unique<atom>();
                //new_atom = new atom;
                new_atom->element = element;
                new_atom->index = atom_index;
                new_atom->coords[0] = xcoord;
                new_atom->coords[1] = ycoord;
                new_atom->coords[2] = zcoord;
                this->atoms.push_back(std::move(new_atom));
                //this->atoms.push_back(new_atom);
                
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


/*
void Structure::get_bond_matrix(){
    int i, j, dim;
    int valence, max_valence;
    int bond_order;
    atom atom_i;
    atom atom_j;

    dim = (this->n_atoms-1)*this->n_atoms/2;
    this->bond_matrix = new int[dim];

    for (i = 0; i < this->n_atoms-1; i++){
        valence = 0;
        atom_i = this->atoms[i];
        max_valence = max_valences[get_element(atom_i.element)];
        for (j = i + 1; j < this->n_atoms; j++){
            if (valence >= max_valence){
                break;
            }
            atom_j = this->atoms[j];
            bond_order = this->get_bond_order(atom_i, atom_j);
            this->bond_matrix[hash_bond_matrix(i, j, this->n_atoms)] = bond_order;
        }
    }

    return;
}
*/


void Structure::get_bonds(){
    int i, j;
    int valence, max_valence;
    int bond_order;
    int atom_i;
    int atom_j;
    bond new_bond;
    std::string element1, element2;
    int terminal_counter;
    
    for (i = 0; i < this->n_atoms-1; i++){
        valence = 0;
        terminal_counter = 0;
        //atom_i = &(this->atoms[i]);
        max_valence = max_valences[this->atoms[i]->element];
        for (j = i + 1; j < this->n_atoms; j++){
            if (valence >= max_valence){
                break;
            }
            //atom_j = &(this->atoms[j]);
            bond_order = this->get_bond_order(i, j);
            if (bond_order){
                this->atoms[i]->bond_partners.push_back(j);
                this->atoms[j]->bond_partners.push_back(i);
                new_bond.atom_index1 = i;
                new_bond.atom_index2 = j;
                new_bond.bond_order = bond_order;
                this->bonds.push_back(new_bond);
                if (is_terminal_atom(this->atoms[j]->element)){
                    terminal_counter++;
                }
            }
        }
        if (terminal_counter == max_valence-1 && !is_terminal_atom(this->atoms[i]->element)){
            this->atoms[i]->core_of_terminal_group = true;
        }
    }

    return;
}


int Structure::get_bond_order(int i, int j)
{
    int bond_order = 0;
    float tolerance = 0.08;
    double single_bond = -1000.0;
    double double_bond = -1000.0;
    double triple_bond = -1000.0;

    float valence_radius_single1 = valence_radii_single[this->atoms[i]->element];
    float valence_radius_single2 = valence_radii_single[this->atoms[j]->element];
    float valence_radius_double1 = valence_radii_double[this->atoms[i]->element];
    float valence_radius_double2 = valence_radii_double[this->atoms[j]->element];
    float valence_radius_triple1 = valence_radii_triple[this->atoms[i]->element];
    float valence_radius_triple2 = valence_radii_triple[this->atoms[j]->element];

    double distance = sqrt(
        pow((this->atoms[i]->coords[0] - this->atoms[j]->coords[0]), 2) +
        pow((this->atoms[i]->coords[1] - this->atoms[j]->coords[1]), 2) +
        pow((this->atoms[i]->coords[2] - this->atoms[j]->coords[2]), 2)
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