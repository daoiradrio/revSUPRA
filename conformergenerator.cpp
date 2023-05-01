#include "conformergenerator.hpp"



ConformerGenerator::ConformerGenerator(){};


ConformerGenerator::~ConformerGenerator(){};


void ConformerGenerator::generate_conformers(Structure* molecule){
    this->get_torsions(&(molecule->bonds));
    std::cout << this->central_torsions.size() << std::endl;
    std::cout << this->terminal_torsions.size() << std::endl;
    std::cout << this->methyl_torsions.size() << std::endl;
    return;
}


void ConformerGenerator::get_torsions_based_on_bond_matrix(){

    return;
}


void ConformerGenerator::get_torsions(std::vector<bond>* bonds){
    atom atom1, atom2;
    std::string element1, element2;
    int valence1, valence2;
    int sum_valences1, sum_valences2;
    int i;

    for (bond bond: *bonds){
        if (bond.bond_order == 1){
            atom1 = bond.atom1;
            atom2 = bond.atom2;
            if (atom1.bond_partners.size() == 1 || atom2.bond_partners.size() == 1){
                continue;
            }
            sum_valences1 = 0;
            for (atom bond_partner: atom1.bond_partners){
                if (max_valences[bond_partner.element] == 1){
                    sum_valences1++;
                }
            }
            sum_valences2 = 0;
            for (atom bond_partner: atom2.bond_partners){
                if (max_valences[bond_partner.element] == 1){
                    sum_valences2++;
                }
            }
            if (sum_valences1 == max_valences[atom1.element]-1){
                if (get_element(atom1.label) == "C"){
                    this->methyl_torsions.push_back(bond);
                }
                else{
                    this->terminal_torsions.push_back(bond);
                }
            }
            else if (sum_valences2 == max_valences[atom2.element]-1){
                if (get_element(atom1.label) == "C"){
                    this->methyl_torsions.push_back(bond);
                }
                else{
                    this->terminal_torsions.push_back(bond);
                }
            }
            else{
                this->central_torsions.push_back(bond);
            }
        }
    }

    return;
}