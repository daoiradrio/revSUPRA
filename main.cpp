#include "structure.cpp"
#include "conformergenerator.cpp"
#include "analyzer.cpp"


int main(){
    std::string filepath = "inputfiles/Propylcyclohexan.xyz";

    Structure mol;
    mol.get_structure(filepath);

    ConformerGenerator gen(mol);
    gen.generate_conformers();

    Analyzer analyzer;
    analyzer.read_xyz(filepath);

    /*
    // CHECK ATOMS AND RESPECTIVE BOND PARTNERS
    std::cout << "______________" << std::endl;
    int atom;
    for (int i = 0; i < gen.mol->atoms.size(); i++){
        std::cout << gen.mol->atoms[i]->element << gen.mol->atoms[i]->index << std::endl;
        for (int j = 0; j < gen.mol->atoms[i]->bond_partners.size(); j++){
            atom = gen.mol->atoms[i]->bond_partners[j];
            std::cout << gen.mol->atoms[atom]->element << gen.mol->atoms[atom]->index << " ";
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