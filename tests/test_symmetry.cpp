#include <gtest/gtest.h>

#include <structure.hpp>
#include <symmetry.hpp>

#include <string>
#include <vector>
#include <memory>
#include <stdlib.h>


TEST(SymmetryTest, DetectRotationalSymmetry)
{   
    int i;

    Symmetry Sym;

    std::shared_ptr<Structure> Alanin_ptr  = std::make_shared<Structure>();
    std::shared_ptr<Structure> Benzol_ptr  = std::make_shared<Structure>();
    std::shared_ptr<Structure> Tyrosin_ptr = std::make_shared<Structure>();

    std::vector<int> alanin_atoms;
    std::vector<int> alanin_methyl_atoms;
    std::vector<int> benzol_atoms;
    std::vector<int> tyrosin_atoms;
    std::vector<int> tyrosin_ring_atoms;

    /*std::string alanin_file  = "/revSUPRA/inputfiles/Alanin.xyz";
    std::string benzol_file  = "/revSUPRA/inputfiles/Benzol.xyz";
    std::string tyrosin_file = "/revSUPRA/inputfiles/Tyrosin.xyz";*/
    std::string home_path = getenv("HOME");
    std::string alanin_file = home_path + "/revSUPRA/inputfiles/Alanin.xyz";
    std::string benzol_file = home_path + "/revSUPRA/inputfiles/Benzol.xyz";
    std::string tyrosin_file = home_path + "/revSUPRA/inputfiles/Tyrosin.xyz";
    

    Alanin_ptr->read_xyz(alanin_file);
    for (i = 0; i < Alanin_ptr->n_atoms; i++){
        alanin_atoms.push_back(i);
    }
    alanin_methyl_atoms = {2, 3, 4, 5};

    EXPECT_EQ(Sym.rot_sym_along_bond(Alanin_ptr, alanin_atoms, 0, 1, 2), false);
    EXPECT_EQ(Sym.rot_sym_along_bond(Alanin_ptr, alanin_atoms, 1, 0, 2), false);
    EXPECT_EQ(Sym.rot_sym_along_bond(Alanin_ptr, alanin_atoms, 0, 1, 3), false);
    EXPECT_EQ(Sym.rot_sym_along_bond(Alanin_ptr, alanin_methyl_atoms, 0, 2, 3), true);
    EXPECT_EQ(Sym.rot_sym_along_bond(Alanin_ptr, alanin_methyl_atoms, 2, 0, 3), true);
    EXPECT_EQ(Sym.rot_sym_along_bond(Alanin_ptr, alanin_methyl_atoms, 0, 2, 2), false);

    Benzol_ptr->read_xyz(benzol_file);
    for (i = 0; i < Benzol_ptr->n_atoms; i++){
        benzol_atoms.push_back(i);
    }

    EXPECT_EQ(Sym.rot_sym_along_bond(Benzol_ptr, benzol_atoms, 0, 8, 2), true);
    EXPECT_EQ(Sym.rot_sym_along_bond(Benzol_ptr, benzol_atoms, 8, 0, 2), true);

    Tyrosin_ptr->read_xyz(tyrosin_file);
    for (i = 0; i < Tyrosin_ptr->n_atoms; i++){
        tyrosin_atoms.push_back(i);
    }
    tyrosin_ring_atoms = {0, 1, 2, 3, 4, 5};

    EXPECT_EQ(Sym.rot_sym_along_bond(Tyrosin_ptr, tyrosin_atoms, 13, 16, 3), false);
    EXPECT_EQ(Sym.rot_sym_along_bond(Tyrosin_ptr, tyrosin_ring_atoms, 4, 12, 2), true);
}