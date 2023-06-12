#include <structure.hpp>
#include <conformergenerator.hpp>
#include <analyzer.hpp>
#include <symmetry.hpp>

#include <iostream>

#include <algorithm> // sort for torsion groups



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

    for (Bond torsion: confgen.torsions){
        confgen.bond_angles.push_back({});
    }

    int 		        i, j;
    int 		        atom1, atom2;
    std::vector<int>	left_atoms;
    std::vector<int>	right_atoms;
    std::vector<int>	status(confgen.mol->n_atoms, 0);
    Symmetry 		    sym;
    std::vector<int>    rot_atoms;

    confgen.input_coords_mat = confgen.mol->coords;

    for (i = 0; i < confgen.mol->coords.rows(); i++){
	    confgen.input_coords.push_back(confgen.mol->coords.row(i));
    }

    for (Bond& torsion: confgen.torsions){
        atom1 = torsion.atom_index1;
        atom2 = torsion.atom_index2;
        left_atoms.clear();
        std::fill(status.begin(), status.end(), 0);
        left_atoms = confgen.torsion_atom_counter(atom1, atom2, status, left_atoms);
        right_atoms.clear();
        std::fill(status.begin(), status.end(), 0);
        right_atoms = confgen.torsion_atom_counter(atom2, atom1, status, right_atoms);
        if (left_atoms.size() <= right_atoms.size()){
            confgen.torsion_atoms.push_back(std::vector<int>());
            i = confgen.torsion_atoms.size() - 1;
            for (int j: left_atoms){
                confgen.torsion_atoms[i].push_back(j);
            }
        }
        else{
            confgen.torsion_atoms.push_back(std::vector<int>());
            i = confgen.torsion_atoms.size() - 1;
            for (int j: right_atoms){
                confgen.torsion_atoms[i].push_back(j);
            }
        }
        //***1
        std::cout << "Bond: " << atom1 << " " << atom2 << std::endl;
        std::fill(status.begin(), status.end(), 0);
        rot_atoms = confgen.get_torsion_group(atom1, atom2, status);
        torsion.rot_sym1 = sym.rot_sym_along_bond(confgen.mol, rot_atoms, atom1, atom2);
        if (torsion.rot_sym1 > 1){
            torsion.rot_sym_atoms1 = rot_atoms;
        }
        std::cout << "Rotationsordnung links: " << torsion.rot_sym1 << std::endl;
        std::fill(status.begin(), status.end(), 0);
        rot_atoms = confgen.get_torsion_group(atom2, atom1, status);
        torsion.rot_sym2 = sym.rot_sym_along_bond(confgen.mol, rot_atoms, atom1, atom2);
        if (torsion.rot_sym2 > 1){
            torsion.rot_sym_atoms2 = rot_atoms;
        }
        std::cout << "Rotationsordnung rechts: " << torsion.rot_sym2 << std::endl;
        std::cout << std::endl;
        //***2
    }

    std::vector<int> container;
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
