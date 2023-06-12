// unbedingt nötig
#include <gtest/gtest.h> // Error weil VSCOde Pfad nicht findet

#include <structure.hpp>

#include <string>
#include <stdlib.h>



//TEST("Test-Klassen Name" bspw. ListInt, Test-Name){
    // normal C++ code for testing...
    //
    // besonders:
    // - EXPECT_EQ bspw. EXPECT_EQ(list.front(), nullptr) << "Error Message"
    // - EXPECT_ANY_THROW
    // - ASSERT_EQ bricht gesamten Test ab wenn false
    // - ...
//}


TEST(TestStructure, ReadInput)
{   
    Structure mol;
    std::string home_path = getenv("HOME");
    std::string file = home_path + "/revSUPRA/inputfiles/D-Alanin.xyz";

    EXPECT_EQ(mol.read_xyz(file), 1);
    EXPECT_EQ(mol.n_atoms, 13);
    EXPECT_EQ(mol.coords.rows(), 13);
    EXPECT_EQ(mol.coords.cols(), 3);
}


TEST(TestStructure, GetConnectivity)
{
    Structure mol;
    std::string home_path = getenv("HOME");
    std::string file = home_path + "/revSUPRA/inputfiles/D-Alanin.xyz";
    mol.get_structure(file);

    EXPECT_EQ(mol.bonds.size(), 12);
    EXPECT_EQ(mol.bonds[4]->bond_order, 2);
    EXPECT_EQ(mol.atoms[2]->core_of_terminal_group, true);
}
