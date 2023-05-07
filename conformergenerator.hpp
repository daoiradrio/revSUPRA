#ifndef CONFORMERGENERATOR_HPP
#define CONFORMERGENERATOR_HPP

#include "helper.hpp"
#include "structure.hpp"

#include <vector>
#include <memory>
#include <string.h>
#include <cmath>
#include <fstream>
#include <queue>
#include <Eigen/Dense>



class ConformerGenerator{

    public:
        ConformerGenerator(Structure input_mol);
        ~ConformerGenerator();

        void generate_conformers();
    
    private:
        std::shared_ptr<Structure> mol;

        std::vector<bond> torsions;
        std::vector<bond> central_torsions;
        std::vector<bond> methylalike_torsions;
        std::vector<bond> terminal_torsions;
        std::vector<std::vector<int>> torsion_atoms;

        std::vector<double> angle_increments; // IN RAD STATT DEGREE??
        std::vector<double> angles;

        std::vector<Eigen::Vector3d> input_coords;

        void get_torsions();

        void find_cycles();
        void cycle_detection(int current, int last, char status[], int ancestors[]);

        void find_peptidebonds();

        //void selection_menu();
        void generation_setup();
        std::vector<int> torsion_atom_counter(int start, int last, int* status, std::vector<int> container);

        int combinations(std::vector<Eigen::Vector3d> new_coords, int index, int counter);

        bool clashes(std::vector<Eigen::Vector3d> coords);
        bool distant_atoms(int atom1, int atom2);

        void write_xyz(std::vector<Eigen::Vector3d> coords, int n = 0);
};


#endif