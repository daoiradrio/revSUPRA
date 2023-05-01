#include "Molecule.hpp"


Molecule::Molecule(){}


Molecule::~Molecule(){}


void Molecule::get_structure(std::string filename)
{
    std::string element, trash;
    double xcoord, ycoord, zcoord;

    std::ifstream infile(filename);

    int iter = 0;
    std::string new_label = "";
    atom new_atom;
    
    while (infile >> element >> xcoord >> ycoord >> zcoord)
    {   
        new_atom.label = element + std::to_string(iter);
        new_atom.coords = {xcoord, ycoord, zcoord};
        new_atom.connectivity = {};
        Molecule::atoms.push_back(new_atom);
        iter++;
    }

    infile.close();
}


void Molecule::get_structure_ORCA(std::string filename)
{
    std::string element, trash;
    double xcoord, ycoord, zcoord;

    std::ifstream infile(filename);
    infile >> Molecule::atom_number;
    getline(infile, trash);
    getline(infile, trash);

    int iter = 0;
    std::string new_label = "";
    atom new_atom;
    while (infile >> element >> xcoord >> ycoord >> zcoord)
    {   
        new_atom.label = element + std::to_string(iter);
        new_atom.coords = {xcoord, ycoord, zcoord};
        new_atom.connectivity = {};
        Molecule::atoms.push_back(new_atom);
        iter++;
    }

    infile.close();

    int iter1 = 0;
    int iter2 = 1;
    int valence = 0;
    int limit = Molecule::atoms.size();
    atom atom1, atom2;
    while (iter1 < limit)
    {
        atom1 = Molecule::atoms[iter1];
        valence = max_valences[get_element(atom1.label)];
        iter2 = iter1 + 1;
        while (iter2 < limit && valence)
        {
            atom2 = Molecule::atoms[iter2];
            if (Molecule::bond_order(&atom1, &atom2))
            {
                Molecule::atoms[iter1].connectivity.push_back(atom2);
                Molecule::atoms[iter2].connectivity.push_back(atom1);
                valence--;
            }
            iter2++;
        }
        iter1++;
    }
}


//ERWEITERN AUF HÖHERE BINDUNGSORDNUNGEN UM ROTIERBARKEIT ZU PRÜFEN
int Molecule::bond_order(atom* atom1, atom* atom2)
{
    int bond_order = 0;
    float tolerance = 0.08;
    float valence_radius1 = valence_radii[get_element(atom1->label)];
    float valence_radius2 = valence_radii[get_element(atom2->label)];

    double single_bond = valence_radius1 + valence_radius2 + tolerance;

    double distance = sqrt(
        pow((atom1->coords[0] - atom2->coords[0]), 2) +
        pow((atom1->coords[1] - atom2->coords[1]), 2) +
        pow((atom1->coords[2] - atom2->coords[2]), 2)
    );

    if (distance <= single_bond)
    {
        bond_order = 1;
    }

    return bond_order;
}


void Molecule::clear()
{
    Molecule::atom_number = 0;
    Molecule::atoms.clear();
}