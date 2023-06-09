#ifndef OPTIMIZER_HPP
#define OPTIMIZER_HPP

#include <structure.hpp>

#include <string>
#include <fstream>
#include <memory>
#include <fstream>
#include <Eigen/Dense>



class Optimizer{

    public:

        int uff_optimization(std::string path, std::string xyz_file, int index = -1);

    private:

};



#endif