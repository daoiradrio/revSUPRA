#include <gtest/gtest.h>

#include <analyzer.hpp>
#include <structure.hpp>

#include <string>
#include <stdlib.h>



TEST(TestAnalyzer, DetectDoubles)
{
    Analyzer analyzer;

    Structure mol1;
    Structure mol2;
    Structure mol3;
    Structure mol4;
    Structure mol5;

    std::string home_path = getenv("HOME");
    std::string file1 = home_path + "/revSUPRA/tests/testcases/Tyrosin.xyz";
    std::string file2 = home_path + "/revSUPRA/tests/testcases/Tyrosin_double.xyz";
    std::string file3 = home_path + "/revSUPRA/tests/testcases/Tyrosin_phenyl_rotated_180_deg.xyz";
    std::string file4 = home_path + "/revSUPRA/tests/testcases/Alanin.xyz";
    std::string file5 = home_path + "/revSUPRA/tests/testcases/D-Alanin.xyz";

    mol1.read_xyz(file1);
    mol2.read_xyz(file2);
    mol3.read_xyz(file3);
    mol4.read_xyz(file4);
    mol5.read_xyz(file5);

    EXPECT_EQ(analyzer.doubles(mol1, mol2), true);
    EXPECT_EQ(analyzer.doubles(mol2, mol1), true);
    EXPECT_EQ(analyzer.doubles(mol1, mol3), true);
    EXPECT_EQ(analyzer.doubles(mol2, mol3), true);
    EXPECT_EQ(analyzer.doubles(mol4, mol5), false);
    EXPECT_EQ(analyzer.doubles(mol5, mol4), false);
}