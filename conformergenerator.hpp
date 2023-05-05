#ifndef CONFORMERGENERATOR_HPP
#define CONFORMERGENERATOR_HPP

#include "helper.hpp"
#include "structure.hpp"

#include <vector>
#include <memory>
#include <string.h>
#include <cmath>



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

        std::vector<int> angle_increments; // IN RAD STATT DEGREE??

        void get_torsions();
        void find_cycles();
        void cycle_detection(int current, int last, char status[], int ancestors[]);
        void find_peptidebonds();
};


#endif