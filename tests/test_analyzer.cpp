#include <gtest/gtest.h>

#include <analyzer.hpp>
#include <structure.hpp>

#include <string>
#include <stdlib.h>



TEST(TestAnalyzer, DetectDoubles)
{
    Analyzer analyzer;

    std::string home_path = getenv("HOME");
    std::string file1 = home_path + "/revSUPRA/tests/testcases/Tyrosin.xyz";
    std::string file2 = home_path + "/revSUPRA/tests/testcases/Tyrosin_double.xyz";
    std::string file3 = home_path + "/revSUPRA/tests/testcases/Tyrosin_phenyl_rotated_180_deg.xyz";
    std::string file4 = home_path + "/revSUPRA/tests/testcases/Alanin.xyz";
    std::string file5 = home_path + "/revSUPRA/tests/testcases/Alanin_methyl_rotated_60_deg.xyz";

    EXPECT_EQ(analyzer.doubles(file1, file2), true);
    EXPECT_EQ(analyzer.doubles(file2, file1), true);
    EXPECT_EQ(analyzer.doubles(file1, file3), true);
    EXPECT_EQ(analyzer.doubles(file2, file3), true);
    EXPECT_EQ(analyzer.doubles(file4, file5), false);
    EXPECT_EQ(analyzer.doubles(file5, file4), false);
}


TEST(TestAnalyzer, DetectDoublesWithoutMethyl)
{
    Analyzer analyzer;

    std::string home_path = getenv("HOME");
    std::string file1 = home_path + "/revSUPRA/tests/testcases/Alanin.xyz";
    std::string file2 = home_path + "/revSUPRA/tests/testcases/Alanin_methyl_rotated_60_deg.xyz";

    EXPECT_EQ(analyzer.doubles(file1, file2, true), true);
}