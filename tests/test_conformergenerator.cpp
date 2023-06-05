#include <gtest/gtest.h>

#include <conformergenerator.hpp>
#include <structure.hpp>

#include <string>
#include <stdlib.h>


TEST(TestConformerGenerator, test)
{
    Structure propylcyclohexan;
    std::string home_path = getenv("HOME");
    std::string ala_ala_file = home_path + "/revSUPRA/inputfiles/Ala-Ala.xyz";
    std::string propylcyclohexan_file = home_path + "/revSUPRA/inputfiles/Propylcyclohexan.xyz";
    propylcyclohexan.get_structure(propylcyclohexan_file);

    ConformerGenerator propylcyclohexan_generator(propylcyclohexan);
    EXPECT_NE(propylcyclohexan_generator.mol, nullptr);

    propylcyclohexan_generator.get_torsions();
    EXPECT_EQ(propylcyclohexan_generator.central_torsions.size(), 8);

    propylcyclohexan_generator.find_cycles();
    EXPECT_EQ(propylcyclohexan_generator.central_torsions.size(), 2);

    //gen.find_peptidebonds();
}