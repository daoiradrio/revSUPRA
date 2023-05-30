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
    int atom_index1;
    int atom_index2;
    int bond_order;
};


/*struct SymmetryElement{
    std::vector<int>    transform;
    int                 order;
    int                 nparam; // NECESSARY??
    double              maxdev;
    double              distance;
    std::vector<double> normal;
    std::vector<double> direction;

    SymmetryElement(): normal(DIMENSION, 0.0), direction(DIMENSION, 0.0) {}
    ~SymmetryElement(){}

    //std::function<void(
    //    std::shared_ptr<SymmetryElement>,
    //    std::shared_ptr<Atom>,
    //    std::shared_ptr<Atom>
    //)> transform_atom;

    void(Symmetry::*transform_atom)(
        std::shared_ptr<SymmetryElement>,
        std::shared_ptr<Atom>,
        std::shared_ptr<Atom>
    );
};*/



#endif
