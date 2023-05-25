#ifndef HELPER_HPP
#define HELPER_HPP

#include <vector>
#include <unordered_map>
#include <cmath>
#include <functional>
#include <memory>

#define DIMENSION        3
#define TolerancePrimary 5e-2
#define ToleranceSame    1e-3



extern std::unordered_map<std::string, int>     element_numbers;
extern std::unordered_map<std::string, double>  valence_radii_single;
extern std::unordered_map<std::string, double>  valence_radii_double;
extern std::unordered_map<std::string, double>  valence_radii_triple;
extern std::unordered_map<std::string, int>     max_valences;


std::string get_element(std::string label);


//int hash_bond_matrix(int row_index, int column_index, int n_atoms){
//    return (2*row_index*n_atoms - row_index - pow(row_index, 2))/2 + column_index - row_index - 1;
//}


bool is_terminal_atom(std::string element);


struct ATOM{
    std::string         element;
    int                 pse_num;
    int                 index;
    std::vector<double> coords;
    std::vector<int>    bond_partners;
    bool                core_of_terminal_group = false;

    ATOM(): coords(DIMENSION, 0.0) {}
    ~ATOM(){}
};


struct BOND{
    int atom_index1;
    int atom_index2;
    int bond_order;
};


struct SYMMETRY_ELEMENT{
    std::vector<int>    transform; // NECESSARY??
    int                 order;
    int                 nparam; // NECESSARY??
    double              maxdev;
    double              distance;
    std::vector<double> normal;
    std::vector<double> direction;

    SYMMETRY_ELEMENT(): normal(DIMENSION, 0.0), direction(DIMENSION, 0.0) {}
    ~SYMMETRY_ELEMENT(){}

    std::function<void(
        std::shared_ptr<struct SYMMETRY_ELEMENT>,
        std::shared_ptr<ATOM> from,
        std::shared_ptr<ATOM> to
    )> transform_atom;
};



#endif