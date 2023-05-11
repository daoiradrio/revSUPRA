#ifndef ANALYZER_HPP
#define ANALYZER_HPP

#include <structure.hpp>

#include <Eigen/Dense>
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

        float rmsd();
        void remove_doubles();
        void kabsch(Eigen::Matrix3d coords1, Eigen::Matrix3d coords2);
        void read_xyz(std::string filepath);

    private:

};


#endif
