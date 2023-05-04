#ifndef ANALYZER_HPP
#define ANALYZER_HPP

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
        std::vector<std::vector<float>> coords;

        float rmsd();
        //void kabsch();
        void read_xyz(std::string filepath);

    private:

};


#endif