#ifndef ROTATIONAXIS_HPP
#define ROTATIONAXIS_HPP

#include <utils.hpp>

#include <memory>
#include <cmath>
#include <Eigen/Dense>



class RotationAxis{

    public:
        Eigen::Vector3d     from;
        Eigen::Vector3d     to;
        Eigen::Vector3d     axis;

        RotationAxis(std::shared_ptr<Atom> from_atom, std::shared_ptr<Atom> to_atom);
	    RotationAxis(std::vector<double> from_coords, std::vector<double> to_coords);
	    RotationAxis(Eigen::Vector3d from_coords, Eigen::Vector3d to_coords);
        ~RotationAxis(){}

        Eigen::Vector3d         rotate_atom(Eigen::Vector3d coords, double deg);
        std::vector<double>     rotate_atom(std::vector<double> coords, double deg);
        std::shared_ptr<Atom>   rotate_atom(std::shared_ptr<Atom> atom, double deg);
        static Eigen::Vector3d  rotate_atom(
                                    const Eigen::Vector3d from_coords,
                                    const Eigen::Vector3d to_coords,
                                    const Eigen::Vector3d coords,
                                    const double          deg
                                );

    private:

};



#endif
