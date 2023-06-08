#ifndef OPTIMIZER_HPP
#define OPTIMIZER_HPP

#include <structure.hpp>

#include <string>
#include <fstream>
#include <memory>
#include <Eigen/Dense>



class Optimizer{

    public:
        std::string workdir_name = "running_opt";

        void uff_optimization(std::string workdir, std::string xyz_file);

    private:

};



#endif