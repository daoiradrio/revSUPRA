#ifndef ROTATIONAXIS_HPP
#define ROTATIONAXIS_HPP

#include <utils.hpp>

#include <memory>
#include <cmath>
#include <Eigen/Dense>



class RotationAxis{

    public:
        int                 order;

        Eigen::Vector3d     from;
        Eigen::Vector3d     to;
        Eigen::Vector3d     axis;

        RotationAxis(std::shared_ptr<Atom> from_atom, std::shared_ptr<Atom> to_atom, int order);
        ~RotationAxis(){}

        Eigen::Vector3d         rotate_atom(Eigen::Vector3d coords);
        std::vector<double>     rotate_atom(std::vector<double> coords);
        std::shared_ptr<Atom>   rotate_atom(std::shared_ptr<Atom> atom);

    private:

};



#endif