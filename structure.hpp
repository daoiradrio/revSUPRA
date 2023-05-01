#include "helper.hpp"



class Structure{

    public:
        Structure();
        ~Structure();

        int n_atoms;
        std::vector<atom> atoms;
        std::vector<bond> bonds;
        int* bond_matrix;

        void get_structure(std::string filepath);
    
    private:
        void read_xyz(std::string filepath);
        void get_bond_matrix();
        int hash_bond_matrix(int row_index, int column_index);
        int get_bond_order(atom* atom1, atom* atom2);
};