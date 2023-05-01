#pragma once

#include <vector>
#include <fstream>
#include <unordered_map>
#include <cmath>

#include "helper.hpp"


class Molecule
{
    public:
        Molecule();
        ~Molecule();

        int atom_number;
        std::vector<atom> atoms;

        void get_structure(std::string filename);
        void get_structure_ORCA(std::string filename);
        
        void clear();
    
    private:
        int bond_order(atom* atom1, atom* atom2);
        //void clear();
};