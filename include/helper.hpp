#ifndef HELPER_HPP
#define HELPER_HPP

#include <vector>
#include <unordered_map>
#include <cmath>



extern std::unordered_map<std::string, int> element_numbers;
extern std::unordered_map<std::string, double> valence_radii_single;
extern std::unordered_map<std::string, double> valence_radii_double;
extern std::unordered_map<std::string, double> valence_radii_triple;
extern std::unordered_map<std::string, int> max_valences;


std::string get_element(std::string label);


//int hash_bond_matrix(int row_index, int column_index, int n_atoms){
//    return (2*row_index*n_atoms - row_index - pow(row_index, 2))/2 + column_index - row_index - 1;
//}


bool is_terminal_atom(std::string element);


struct atom{
    std::string element;
    int index;
    std::vector<double> coords;
    std::vector<int> bond_partners;
    bool core_of_terminal_group = false;
};


struct bond{
    int atom_index1;
    int atom_index2;
    int bond_order;
};


#endif