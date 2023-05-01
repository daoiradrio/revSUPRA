//#include "structure.cpp"
#include "conformergenerator.cpp"



int main(){
    std::string filepath = ".xyz-Inputdateien/Alanin.xyz";

    Structure mol;
    ConformerGenerator gen;

    mol.get_structure(filepath);
    //gen.generate_conformers(&mol);

    /*std::cout << "Central torsions" << std::endl;
    for (bond b: gen.central_torsions){
        std::cout << b.atom1.label << "  " << b.atom2.label << std::endl;
    }
    std::cout << "Terminal torsions" << std::endl;
    for (bond b: gen.terminal_torsions){
        std::cout << b.atom1.label << "  " << b.atom2.label << std::endl;
    }
    std::cout << "Methyl torsions" << std::endl;
    for (bond b: gen.methyl_torsions){
        std::cout << b.atom1.label << "  " << b.atom2.label << std::endl;
    }*/

    return 0;
}