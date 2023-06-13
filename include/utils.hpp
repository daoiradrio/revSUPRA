#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <unordered_map>
#include <cmath>
#include <functional>
#include <memory>

#define DIMENSION           3
#define TolerancePrimary    5e-2
#define ToleranceSame       1e-3
#define ToleranceFinal      1e-4
#define MaxOptCycles        200
#define MaxOptStep          5e-1
#define MinOptStep          1e-7
#define OptChangeThreshold  1e-10
#define GradientStep        1e-7
#define OptChangeHits       5
#define BUFFER_SIZE         256
#define SUCCESS_EXIT        1
#define FAIL_EXIT           0



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


struct Atom{
    std::string         element;
    int                 pse_num;
    int                 index;
    std::vector<double> coords;
    std::vector<int>    bond_partners;
    bool                core_of_terminal_group = false;

    Atom(): coords(DIMENSION, 0.0) {}
    ~Atom(){}
};


struct Bond{
    std::shared_ptr<Atom>   atom1;
    std::shared_ptr<Atom>   atom2;
    int                     bond_order = 0;
    int                     rot_sym1   = 1;
    int                     rot_sym2   = 1;
    std::vector<int>        rot_sym_atoms1;
    std::vector<int>        rot_sym_atoms2;
};


struct Torsion{
    
    std::shared_ptr<Bond>   bond;
    std::vector<int>        rot_atoms;
};



#endif
