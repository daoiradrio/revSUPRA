#ifndef SYMMETRYELEMENT_HPP
#define SYMMETRYELEMENT_HPP

#include <utils.hpp>

#include <iostream>
#include <functional>
#include <memory>
#include <vector>



/* WAS IST DER GOOD PRACTICE WEG, UM KLASSEN MITEINANDER INTERAGIEREN ZU LASSEN? */



class SymmetryElement{
    public:
        std::vector<int>    transform;
        int                 order;
        double              maxdev;
        double              distance;
        std::vector<double> normal;
        std::vector<double> direction;

        SymmetryElement(): normal(DIMENSION, 0.0), direction(DIMENSION, 0.0) {}
        ~SymmetryElement(){}

        virtual int transform_atom(std::shared_ptr<Atom> from, std::shared_ptr<Atom> to);

    private:

};



#endif