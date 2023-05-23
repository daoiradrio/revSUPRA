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
#include <limits>
#include <algorithm>
#include <cmath>
#include <utility>


class Analyzer{

    public:
        Analyzer();
        ~Analyzer();

        std::vector<std::string> elements;
        Eigen::MatrixX3d coords;

        double rmsd(Eigen::MatrixX3d coords1, Eigen::MatrixX3d coords2);
        void remove_doubles(std::string filepath, std::string filename);
        void extract_energies(std::string folderpath, std::string foldername);
        void divide_and_conquer_remove_doubles(std::string filepath, std::string filename);
        bool doubles(Structure& struc1, Structure& struc2);

        std::vector<std::pair<double, int>> container;
        std::vector<double> energies;

    private:
        //std::vector<std::pair<double, int>> container;

        static bool sort_func(std::pair<double, int> a, std::pair<double, int> b){return (a.first < b.first);};

};


#endif
