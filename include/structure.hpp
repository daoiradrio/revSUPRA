#ifndef STRUCTURE_HPP
#define STRUCTURE_HPP

#include <utils.hpp>

#include <vector>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <Eigen/Dense>



class Structure{

    public:
        Structure();
        ~Structure();

        int                                 n_atoms;
        std::vector<std::shared_ptr<Atom>>  atoms;
        std::vector<Bond>  bonds;
        Eigen::MatrixX3d                    coords;

        void    get_structure(std::string filepath);
        //void read_xyz(std::string filepath);
        int     read_xyz(std::string filepath);
        //void read_xyz(std::string filepath);
        //void get_bonds();
    
    private:
        void    get_connectivity();
        int     get_bond_order(int i, int j);

        //int* bond_matrix;
        //void get_bond_matrix();
};


#endif
