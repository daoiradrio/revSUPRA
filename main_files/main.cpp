#include <structure.hpp>
#include <conformergenerator.hpp>
#include <analyzer.hpp>

#include <iostream>


int main(int argc, char **argv){
    std::string filename;
    /*if (argc == 2){
        filename = argv[1];
    }
    else{
        std::cout << "Type in .xyz file name: ";
        std::cin >> filename;
    }
    std::string filepath = "inputfiles/" + filename;*/

    std::string filepath = "/home/baum/SUPRA/inputfiles/Tyrosin.xyz";

    Structure mol;
    mol.get_structure(filepath);

    ConformerGenerator gen(mol);
    gen.generate_conformers();

    /*
    std::string clash_filepath = "debug_files/clash_structure1.xyz";
    Structure clash_mol;
    clash_mol.get_structure(clash_filepath);
    std::vector<Eigen::Vector3d> coords;
    for (int i = 0; i < clash_mol.n_atoms; i++){
        Eigen::Vector3d coord(clash_mol.atoms[i]->coords.data());
        coords.push_back(coord);
    }
    if (gen.clashes(coords)){
        std::cout << "CLASHES!" << std::endl;
    }
    else{
        std::cout << "NO CLASHES." << std::endl;
    }
    */

    //Analyzer analyzer;
    //analyzer.read_xyz(filepath);

    /*
    // CHECK ATOMS AND RESPECTIVE BOND PARTNERS
    std::cout << "______________" << std::endl;
    int atom;
    for (int i = 0; i < gen->mol->atoms.size(); i++){
        std::cout << gen->mol->atoms[i]->element << gen->mol->atoms[i]->index << std::endl;
        for (int j = 0; j < gen->mol->atoms[i]->bond_partners.size(); j++){
            atom = gen->mol->atoms[i]->bond_partners[j];
            std::cout << gen->mol->atoms[atom]->element << gen->mol->atoms[atom]->index << " ";
        }
        std::cout << std::endl;
    }
    */

    /*
    // CHECK NUMBER OF DIFFERENT TYPES OF TORSIONS
    std::cout << gen.central_torsions.size() << std::endl;
    std::cout << gen.terminal_torsions.size() << std::endl;
    std::cout << gen.methylalike_torsions.size() << std::endl;
    */

    return 0;
}