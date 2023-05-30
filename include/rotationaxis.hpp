#ifndef ROTATIONAXIS_HPP
#define ROTATIONAXIS_HPP

#include <utils.hpp>
#include <symmetryelement.hpp>

#include <memory>
#include <cmath>



class RotationAxis: public SymmetryElement{

    public:
        using SymmetryElement::transform_atom;
        int transform_atom(std::shared_ptr<Atom> from, std::shared_ptr<Atom> to);

    private:

};



#endif