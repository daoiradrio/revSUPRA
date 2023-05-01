#include <iostream>
#include "Molecule.cpp"

#include "helper.hpp"
#include <cmath>
#include <stdlib.h>

#include <Eigen/Dense>


int main()
{   
    Molecule molecule1;
    Molecule molecule2;

    std::string folder = "AlaninA60/";
    std::string file = "conformer";
    std::string extension = ".xyz";

    int limit = 1295;
    int counter = 0;
    for (int i = 1; i <= limit - 1; i++)
    {
        molecule1.clear();
        molecule1.get_structure_ORCA(folder + file + std::to_string(i) + extension);

        if (molecule1.atoms.size() == 0)
        {
            continue;
        }  

        for (int j = i + 1; j <= limit; j++)
        {
            molecule2.clear();
            molecule2.get_structure_ORCA(folder + file + std::to_string(j) + extension);

            if (molecule2.atoms.size() == 0)
            {
                continue;
            }
            
            if (quaternion_rmsd(molecule1.atoms, molecule2.atoms) <= 0.1)
            {
                counter++;
                break;
            }
        }
    }

    std::cout << limit + 1 - counter << " individuelle Konformere" << std::endl;

    return 0;
}