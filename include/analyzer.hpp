#ifndef ANALYZER_HPP
#define ANALYZER_HPP

#include <structure.hpp>
#include <hungarian.hpp>

#include <Eigen/Dense>
#include <Eigen/SVD>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>


class Analyzer{

    public:
        Analyzer();
        ~Analyzer();

        std::vector<std::string> elements;
        Eigen::MatrixX3d coords;

        double rmsd(Eigen::MatrixX3d coords1, Eigen::MatrixX3d coords2);
        void remove_doubles(std::string filepath, std::string filename, int n_files);

    private:

};


#endif
