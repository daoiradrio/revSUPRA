//#include "structure.cpp"
#include "conformergenerator.cpp"



int main(){
    std::string filepath = "inputfiles/Alanin.xyz";

    Structure mol;
    std::unique_ptr<Structure> ptr(&mol);
    ConformerGenerator gen(&ptr);

    mol.get_structure(filepath);
    //gen.generate_conformers();

    /*for (bond* torsion: gen.central_torsions){
        std::cout << torsion->atom1->label << " " << torsion->atom2->label << std::endl;
    }

    for (bond* rb: gen.ring_bonds){
        std::cout << rb->atom1->label << " " << rb->atom2->label << std::endl;
    }*/

    return 0;
}