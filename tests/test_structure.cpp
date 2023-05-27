// unbedingt nötig
#include <gtest/gtest.h> // Error weil VSCOde Pfad nicht findet

//#include <structure.hpp>

//#include <string>



//TEST("Test-Klassen Name" bspw. ListInt, Test-Name){
    // normal C++ code for testing...
    //
    // besonders:
    // - EXPECT_EQ bspw. EXPECT_EQ(list.front(), nullptr) << "Error Message"
    // - EXPECT_ANY_THROW
    // - ASSERT_EQ bricht gesamten Test ab wenn false
    // - ...
//}


TEST(StructureTest, ReadInput)
{   
    //Structure mol;
    //std::string file = "Alanin.xyz";
    //EXPECT_EQ(mol.read_xyz(file), 1);
    EXPECT_EQ(7*6, 42); // PLACEHOLDER, DELETE LATER
}
