#include <gtest/gtest.h>

#include <conformergenerator.hpp>
#include <structure.hpp>

#include <string>


TEST(TestConformerGenerator, test)
{
    Structure propylcyclohexan;
    std::string ala_ala_file = "/home/baum/revSUPRA/inputfiles/Ala-Ala.xyz";
    std::string propylcyclohexan_file = "/home/baum/revSUPRA/inputfiles/Propylcyclohexan.xyz";
    propylcyclohexan.get_structure(propylcyclohexan_file);

    ConformerGenerator propylcyclohexan_generator(propylcyclohexan);
    EXPECT_NE(propylcyclohexan_generator.mol, nullptr);

    propylcyclohexan_generator.get_torsions();
    EXPECT_EQ(propylcyclohexan_generator.central_torsions.size(), 8);

    propylcyclohexan_generator.find_cycles();
    EXPECT_EQ(propylcyclohexan_generator.central_torsions.size(), 2);

    //gen.find_peptidebonds();
}