#include <structure.hpp>
#include <conformergenerator.hpp>
#include <analyzer.hpp>
#include <symmetry.hpp>

#include <iostream>

#include <algorithm> // min, max



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
    //confgen.generate_conformers();
    confgen.get_torsions();
    confgen.find_cycles();
    confgen.find_peptidebonds();
    confgen.selection_menu();
    confgen.generation_setup();

    for (const std::shared_ptr<Torsion>& torsion: confgen.torsions){
        confgen.bond_angles.push_back({{}, {}});
    }

    Symmetry            sym;
    int                 i, j, k;
    int                 increment_iter;
    int                 order_left;
    int                 order_right;
    std::vector<int>    left_torsion_group;
    std::vector<int>    right_torsion_group;
    std::vector<int>	status(confgen.mol->n_atoms, 0);
    std::vector<int>    torsion_done(confgen.torsions.size(), 0);

    for (std::shared_ptr<Torsion>& torsion: confgen.torsions){
        std::fill(status.begin(), status.end(), 0);
        left_torsion_group = confgen.get_torsion_group(
            torsion->bond->atom1->index,
            torsion->bond->atom2->index,
            status
        );
        torsion->bond->rot_sym1 = sym.rot_sym_along_bond(
            confgen.mol,
            left_torsion_group,
            torsion->bond->atom1->index,
            torsion->bond->atom2->index
        );
        if (torsion->bond->rot_sym1 > 1){
            torsion->bond->rot_sym_atoms1 = left_torsion_group;
        }
        std::fill(status.begin(), status.end(), 0);
        right_torsion_group = confgen.get_torsion_group(
            torsion->bond->atom2->index,
            torsion->bond->atom1->index,
            status
        );
        torsion->bond->rot_sym2 = sym.rot_sym_along_bond(
            confgen.mol,
            right_torsion_group,
            torsion->bond->atom2->index,
            torsion->bond->atom1->index
        );
        if (torsion->bond->rot_sym2 > 1){
            torsion->bond->rot_sym_atoms2 = right_torsion_group;
        }
    }

    i = 0;
    for (std::shared_ptr<Torsion>& torsion1: confgen.torsions){
        order_left = sym.rot_sym_along_bond(
                        confgen.mol,
                        torsion1->left_atoms,
                        torsion1->bond->atom1->index,
                        torsion1->bond->atom2->index
                     );
        if (order_left > 1){
            torsion1->bond->rot_sym1 = order_left;
            torsion1->bond->rot_sym_atoms1 = torsion1->left_atoms;
            std::cout << "C" << order_left << "-Symmetrie auf Seite von " 
                      << torsion1->bond->atom1->element << torsion1->bond->atom1->index << std::endl;
            increment_iter = 0;
            for (int increment: confgen.angle_increments){
                if (std::max(360/order_left, increment)%std::min(360/order_left, increment) == 0){
                    for (j = 0; j*increment < 360/order_left; j++){
                        confgen.bond_angles[i][increment_iter].push_back(j*increment);
                    }
                }
                else{
                    for (j = 0; j < 360/increment; j++){
                        confgen.bond_angles[i][increment_iter].push_back(j*increment);
                    }
                }
                increment_iter++;
            }
            torsion_done[i] = 1;
            i++;
            continue;
        }
        order_right = sym.rot_sym_along_bond(
                        confgen.mol,
                        torsion1->right_atoms,
                        torsion1->bond->atom1->index,
                        torsion1->bond->atom2->index
                     );
        if (order_right > 1){
            torsion1->bond->rot_sym2 = order_right;
            torsion1->bond->rot_sym_atoms2 = torsion1->right_atoms;
            std::cout << "C" << order_right << "-Symmetrie auf Seite von " 
                      << torsion1->bond->atom2->element << torsion1->bond->atom2->index << std::endl;
            increment_iter = 0;
            for (int increment: confgen.angle_increments){
                if (std::max(360/order_right, increment)%std::min(360/order_right, increment) == 0){
                    for (j = 0; j*increment < 360/order_right; j++){
                        confgen.bond_angles[i][increment_iter].push_back(j*increment);
                    }
                }
                else{
                    for (j = 0; j < 360/increment; j++){
                        confgen.bond_angles[i][increment_iter].push_back(j*increment);
                    }
                }
                increment_iter++;
            }
            torsion_done[i] = 1;
            i++;
            continue;
        }
        j = 0;
        for (std::shared_ptr<Torsion>& torsion2: confgen.torsions){
            //if ((torsion1->bond->atom1->index == torsion2->bond->atom1->index  &&
            //     torsion1->bond->atom2->index == torsion2->bond->atom2->index) ||
            //    (torsion1->bond->atom1->index == torsion2->bond->atom2->index  &&
            //     torsion1->bond->atom2->index == torsion2->bond->atom1->index)
            //){
            if (i == j){
                j++;
                continue;
            }
            if (torsion1->bond->rot_sym_atoms1.size() != 0 || torsion1->bond->rot_sym_atoms1 == torsion2->bond->rot_sym_atoms2){
                if(torsion1->bond->rot_sym_atoms1 == torsion2->bond->rot_sym_atoms1){
                    std::cout << torsion1->bond->atom1->index << " " << torsion1->bond->atom2->index << " und " 
                              << torsion2->bond->atom1->index << " " << torsion2->bond->atom2->index << std::endl;
                    increment_iter = 0;
                    for (int increment: confgen.angle_increments){
                        if (std::max(360/torsion1->bond->rot_sym1, increment)%std::min(360/torsion1->bond->rot_sym1, increment) == 0){
                            for (k = 0; k*increment < 360/torsion1->bond->rot_sym1; k++){
                                confgen.bond_angles[i][increment_iter].push_back(k*increment);
                            }
                        }
                        else{
                            for (k = 0; k < 360/increment; k++){
                                confgen.bond_angles[i][increment_iter].push_back(k*increment);
                            }
                        }
                        increment_iter++;
                    }
                }
                //else if(torsion1->bond->rot_sym_atoms1 == torsion2->bond->rot_sym_atoms2){
                //    std::cout << torsion1->bond->atom1->index << " " << torsion1->bond->atom2->index << " und " 
                //              << torsion2->bond->atom1->index << " " << torsion2->bond->atom2->index << std::endl;
                //}
            }
            if (torsion1->bond->rot_sym_atoms2.size() != 0){
                if(torsion1->bond->rot_sym_atoms2 == torsion2->bond->rot_sym_atoms1){
                    std::cout << torsion1->bond->atom1->index << " " << torsion1->bond->atom2->index << " und " 
                              << torsion2->bond->atom1->index << " " << torsion2->bond->atom2->index << std::endl;
                }
                else if(torsion1->bond->rot_sym_atoms2 == torsion2->bond->rot_sym_atoms2){
                    std::cout << torsion1->bond->atom1->index << " " << torsion1->bond->atom2->index << " und " 
                              << torsion2->bond->atom1->index << " " << torsion2->bond->atom2->index << std::endl;
                }
            }
            j++;
        }
        i++;
        std::cout << std::endl;
    }

    i = 0;
    for (std::shared_ptr<Torsion> torsion: confgen.torsions){
        std::cout << torsion->bond->atom1->index << " " << torsion->bond->atom2->index << std::endl;
        for (j = 0; j < confgen.angle_increments.size(); j++){
            std::cout << confgen.angle_increments[j] << std::endl;
            for (int angle: confgen.bond_angles[i][j]){
                std::cout << angle << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    /*std::vector<int> container;
    for (Bond& torsion1: confgen.torsions){
        std::cout << std::endl;
        container.clear();
        if (torsion1.rot_sym1 == 1 && torsion1.rot_sym2 == 1){
            for (int increment: confgen.angle_increments){
                for (i = 0; i < 360/increment; i++){
                    container.push_back(i*increment);
                }
            }
            confgen.bond_angles.push_back(container);
        }
        else{
            for (Bond& torsion2: confgen.torsions){
                if (torsion1.atom_index1 == torsion2.atom_index1 && torsion1.atom_index2 == torsion2.atom_index2){
                    continue;
                }
                if (torsion1.rot_sym_atoms1.size() != 0){
                    if (torsion1.rot_sym_atoms1 == torsion2.rot_sym_atoms1){
                        std::cout << torsion1.atom_index1 << " " << torsion1.atom_index2 << " und " << torsion2.atom_index1 << " " << torsion2.atom_index2 << std::endl;
                    }
                    if (torsion1.rot_sym_atoms1 == torsion2.rot_sym_atoms2){
                        std::cout << torsion1.atom_index1 << " " << torsion1.atom_index2 << " und " << torsion2.atom_index1 << " " << torsion2.atom_index2 << std::endl;
                    }
                }
                if (torsion1.rot_sym_atoms2.size() != 0){
                    if (torsion1.rot_sym_atoms2 == torsion2.rot_sym_atoms1){
                        std::cout << torsion1.atom_index1 << " " << torsion1.atom_index2 << " und " << torsion2.atom_index1 << " " << torsion2.atom_index2 << std::endl;
                    }
                    if (torsion1.rot_sym_atoms2 == torsion2.rot_sym_atoms2){
                        std::cout << torsion1.atom_index1 << " " << torsion1.atom_index2 << " und " << torsion2.atom_index1 << " " << torsion2.atom_index2 << std::endl;
                    }
                }
            }
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
