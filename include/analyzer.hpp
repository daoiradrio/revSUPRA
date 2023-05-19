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


class Analyzer{

    public:
        Analyzer();
        ~Analyzer();

        std::vector<std::string> elements;
        Eigen::MatrixX3d coords;

        double rmsd(Eigen::MatrixX3d coords1, Eigen::MatrixX3d coords2);
        void remove_doubles(std::string filepath, std::string filename, int n_files);
        void extract_energies(std::string folderpath, std::string foldername, int n_folders);
        void divide_and_conquer_remove_doubles(std::string filepath, std::string filename, int n_files);

        std::vector<std::vector<double>> container;

    private:
        //std::vector<std::vector<double>> container;

        static bool sort_func(std::vector<double> a, std::vector<double> b){return (a[0] < b[0]);};

};


#endif
