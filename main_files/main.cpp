#include <structure.hpp>
#include <conformergenerator.hpp>
#include <analyzer.hpp>
#include <symmetry.hpp>

#include <iostream>



int main(int argc, char **argv){
    std::string filename;
    if (argc == 2){
        filename = argv[1];
    }
    else{
        std::cout << "Type in .xyz file name: ";
        std::cin >> filename;
    }
    std::string filepath = "inputfiles/" + filename;

    Structure mol;
    mol.get_structure(filepath);

    ConformerGenerator confgen(mol);
    confgen.generate_conformers();
    /*confgen.get_torsions();
    confgen.find_cycles();
    confgen.find_peptidebonds();
    confgen.selection_menu();
    confgen.generation_setup();*/

    // PRINT CONNECTIVITY
    /*for (auto atom1: mol.atoms){
        std::cout << atom1->element << atom1->index << ": ";
        for (int atom2: atom1->bond_partners){
            std::cout << mol.atoms[atom2]->element << mol.atoms[atom2]->index << " ";
        }
        std::cout << std::endl;
    }*/

    /*Symmetry                    sym;
    int                         i, j, k;
    std::vector<int>            left_torsion_group;
    std::vector<int>            right_torsion_group;
    std::vector<int>	        status(confgen.mol->n_atoms, 0);
    std::vector<int>            torsion_done(confgen.torsions.size(), 0);
    std::shared_ptr<Torsion>    torsion1;
    std::shared_ptr<Torsion>    torsion2;

    // NECESSARY: algorithm HEADER FOR min, max

    for (std::shared_ptr<Torsion>& torsion: confgen.torsions){
        std::fill(status.begin(), status.end(), 0);
        left_torsion_group = confgen.get_torsion_group(
            torsion->bond->atom1->index,
            torsion->bond->atom2->index,
            status
        );
        torsion->rot_sym1 = sym.rot_sym_along_bond(
            confgen.mol,
            left_torsion_group,
            torsion->bond->atom1->index,
            torsion->bond->atom2->index
        );
        if (torsion->rot_sym1 > 1){
            torsion->rot_sym_atoms1 = left_torsion_group;
        }
        std::fill(status.begin(), status.end(), 0);
        right_torsion_group = confgen.get_torsion_group(
            torsion->bond->atom2->index,
            torsion->bond->atom1->index,
            status
        );
        torsion->rot_sym2 = sym.rot_sym_along_bond(
            confgen.mol,
            right_torsion_group,
            torsion->bond->atom2->index,
            torsion->bond->atom1->index
        );
        if (torsion->rot_sym2 > 1){
            torsion->rot_sym_atoms2 = right_torsion_group;
        }
    }

    for (int increment: confgen.angle_increments){
        confgen.rot_angles.clear();
        for (const std::shared_ptr<Torsion>& torsion: confgen.torsions){
            confgen.rot_angles.push_back({});
        }
        std::fill(torsion_done.begin(), torsion_done.end(), 0);
        for (i = 0; i < confgen.torsions.size(); i++){
            if (torsion_done[i]){
                continue;
            }
            torsion1 = confgen.torsions[i];
            if (torsion1->rot_sym1 == 1 || std::max(360/torsion1->rot_sym1, increment) % std::min(360/torsion1->rot_sym1, increment) != 0){
                if (torsion1->rot_sym2 == 1 || std::max(360/torsion1->rot_sym2, increment) % std::min(360/torsion1->rot_sym2, increment) != 0){
                    for (k = 0; k < 360/increment; k++){
                        confgen.rot_angles[i].push_back(k*increment);
                    }
                    torsion_done[i] = 1;
                    continue;
                }
            }
            if (torsion1->rot_sym_atoms1 == torsion1->rot_atoms1){
                if (torsion1->rot_sym1 > 1){
                    if (std::max(360/torsion1->rot_sym1, increment)%std::min(360/torsion1->rot_sym1, increment) == 0){
                        for (j = 0; j*increment < 360/torsion1->rot_sym1; j++){
                            confgen.rot_angles[i].push_back(j*increment);
                        }
                        torsion_done[i] = 1;
                        continue;
                    }
                }
            }
            if (torsion1->rot_sym_atoms2 == torsion1->rot_atoms2){
                if (torsion1->rot_sym2 > 1){
                    if (std::max(360/torsion1->rot_sym2, increment)%std::min(360/torsion1->rot_sym2, increment) == 0){
                        for (j = 0; j*increment < 360/torsion1->rot_sym2; j++){
                            confgen.rot_angles[i].push_back(j*increment);
                        }
                        torsion_done[i] = 1;
                        continue;
                    }
                }
            }
            for (j = i+1; j < confgen.torsions.size(); j++){
                torsion2 = confgen.torsions[j];
                if (torsion1->rot_sym1 > 1){
                    if (std::max(360/torsion1->rot_sym1, increment)%std::min(360/torsion1->rot_sym1, increment) == 0){
                        if (torsion1->rot_sym_atoms1 == torsion2->rot_sym_atoms1 || torsion1->rot_sym_atoms1 == torsion2->rot_sym_atoms2){
                            if (!torsion_done[j]){
                                for (k = 0; k*increment < 360/torsion1->rot_sym1; k++){
                                    confgen.rot_angles[j].push_back(k*increment);
                                }
                                torsion_done[j] = 1;
                                for (k = 0; k < 360/increment; k++){
                                    confgen.rot_angles[i].push_back(k*increment);
                                }
                                torsion_done[i] = 1;
                                continue;
                            }
                        }
                    }
                }
                if (torsion1->rot_sym2 > 1){
                    if (std::max(360/torsion1->rot_sym2, increment)%std::min(360/torsion1->rot_sym2, increment) == 0){
                        if (torsion1->rot_sym_atoms2 == torsion2->rot_sym_atoms1 || torsion1->rot_sym_atoms2 == torsion2->rot_sym_atoms2){
                            if (!torsion_done[j]){
                                for (k = 0; k*increment < 360/torsion1->rot_sym2; k++){
                                    confgen.rot_angles[j].push_back(k*increment);
                                }
                                torsion_done[j] = 1;
                                for (k = 0; k < 360/increment; k++){
                                    confgen.rot_angles[i].push_back(k*increment);
                                }
                                torsion_done[i] = 1;
                                continue;
                            }
                        }
                    }
                }
            }
        }
        std::cout << "Inkrement " << increment << std::endl;
        std::cout << std::endl;
        i = 0;
        for (std::shared_ptr<Torsion> torsion: confgen.torsions){
            std::cout << torsion->bond->atom1->element << torsion->bond->atom1->index << " "
                      << torsion->bond->atom2->element << torsion->bond->atom2->index << std::endl;
            for (int angle: confgen.rot_angles[i]){
                std::cout << angle << " ";
            }
            std::cout << "\n\n";
            i++;
        }
    }*/

    /*
    vector< vector<double> > costMatrix = {{ 50, 1, 51, 52},
                                           { 1, 50, 51, 52},
                                           { 50, 51, 1, 52},
                                           { 50, 51, 52, 1}};

    HungarianAlgorithm HungAlgo;
    vector<int> assignment;

    double cost = HungAlgo.Solve(costMatrix, assignment);

    for (unsigned int x = 0; x < costMatrix.size(); x++)
        std::cout << x << "," << assignment[x] << "\t";

    std::cout << "\ncost: " << cost << std::endl;
    */

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
