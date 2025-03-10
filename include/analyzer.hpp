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
//#include <limits>
//#include <algorithm>
#include <cmath>
//#include <utility>



class Analyzer{

    public:
        Analyzer();
        ~Analyzer();

        void    remove_doubles(
                    std::string filepath,
                    std::string filename,
                    bool ignore_methyl = false,
                    double rmsd_threshold = 0.1);
        //void    remove_doubles(std::string filepath, double rmsd_threshold = 0.1);
        bool    doubles(
                    std::string file1,
                    std::string file2,
                    bool ignore_methyl = false,
                    double rmsd_threshold = 0.1
                );
        bool    doubles(
                    const Structure& mol1,
                    const Structure& mol2,
                    bool ignore_methyl = false,
                    double rmsd_threshold = 0.1
                );
        void    match_coords(
                    std::vector<std::shared_ptr<Atom>> atoms1,
                    std::vector<std::shared_ptr<Atom>> atoms2,
                    Eigen::MatrixX3d coords1,
                    Eigen::MatrixX3d& coords2
                );
        //void    match_coords_tight();
        void    kabsch(Eigen::MatrixX3d& coords1, Eigen::MatrixX3d& coords2);
        double  rmsd(Eigen::MatrixX3d coords1, Eigen::MatrixX3d coords2);

        //std::vector<std::string> elements;
        //Eigen::MatrixX3d coords;

        //void extract_energies(std::string folderpath, std::string foldername);
        //void divide_and_conquer_remove_doubles(std::string filepath, std::string filename);
        //std::vector<std::pair<double, int>> container;
        //std::vector<double> energies;

    private:
        void remove_methly_atoms(
                const Structure& mol,
                Eigen::MatrixX3d& coords,
                std::vector<std::shared_ptr<Atom>>& atoms
             );

        //std::vector<std::pair<double, int>> container;
        //static bool sort_func(std::pair<double, int> a, std::pair<double, int> b){return (a.first < b.first);};

};



#endif
