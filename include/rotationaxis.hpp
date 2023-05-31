#ifndef ROTATIONAXIS_HPP
#define ROTATIONAXIS_HPP

#include <utils.hpp>
#include <symmetryelement.hpp>

#include <memory>
#include <cmath>
#include <Eigen/Dense>



class RotationAxis: public SymmetryElement{

    public:
        std::shared_ptr<Atom>   axis_from;
        std::shared_ptr<Atom>   axis_to;

        RotationAxis(){this->nparam = 7;}
        ~RotationAxis(){}

        using SymmetryElement::transform_atom;
        std::shared_ptr<Atom> transform_atom(std::shared_ptr<Atom> from);

    private:

};



#endif