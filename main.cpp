//#include "structure.cpp"
#include "conformergenerator.cpp"



int main(){
    std::string filepath = ".xyz-Inputdateien/Alanin.xyz";

    Structure mol;
    ConformerGenerator gen;

    mol.get_structure(filepath);

    for (int i = 0; i < mol.bonds.size(); i++){
        std::cout << mol.bonds[i].atom1.label << " " << mol.bonds[i].atom2.label << std::endl;
        std::cout << mol.bonds[i].atom1.bond_partners.size() << " " << mol.bonds[i].atom2.bond_partners.size() << std::endl;
    }

    gen.generate_conformers(&mol);

    /*for (bond* torsion: gen.central_torsions){
        std::cout << torsion->atom1->label << " " << torsion->atom2->label << std::endl;
    }

    for (bond* rb: gen.ring_bonds){
        std::cout << rb->atom1->label << " " << rb->atom2->label << std::endl;
    }*/

    return 0;
}