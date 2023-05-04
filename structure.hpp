#include "helper.hpp"



class Structure{

    public:
        Structure();
        ~Structure();

        int n_atoms;
        std::vector<std::unique_ptr<atom>> atoms;
        //std::vector<atom*> atoms;
        std::vector<bond> bonds;

        void get_structure(std::string filepath);
    
    private:
        void read_xyz(std::string filepath);
        void get_bonds();
        int get_bond_order(int i, int j);

        //int* bond_matrix;
        //void get_bond_matrix();
};