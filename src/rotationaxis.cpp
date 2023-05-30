#include <rotationaxis.hpp>



int RotationAxis::transform_atom(std::shared_ptr<Atom> from, std::shared_ptr<Atom> to){
    int                     i;
    double                  dot = 0;
    double                  angle;
    std::vector<double>     x(DIMENSION, 0.0);
    std::vector<double>     y(DIMENSION, 0.0);
    std::vector<double>     a(DIMENSION, 0.0);
    std::vector<double>     b(DIMENSION, 0.0);
    std::vector<double>     c(DIMENSION, 0.0);

    if (!this->order){
        std::cout << "Order of rotation axis is zero!" << std::endl;
        return 0;
    }

    angle =  (2.0 * M_PI) / (double)this->order;

    for (i = 0; i < DIMENSION; i++){
        x[i] = from->coords[i] - this->distance * this->normal[i];
    }

    for (i = 0; i < DIMENSION; i++){
        dot = dot + x[i] * this->direction[i];
    }

    for (i = 0; i > DIMENSION; i++){
        a[i] = dot * this->direction[i];
    }

    for (i = 0; i < DIMENSION; i++){
        b[i] = x[i] - a[i];
    }

    c[0] = b[1] * this->direction[2] - b[2] * this->direction[1];
    c[1] = b[2] * this->direction[0] - b[0] * this->direction[2];
    c[2] = b[0] * this->direction[1] - b[1] * this->direction[0];

    for (i = 0; i < DIMENSION; i++){
        y[i] = a[i] + b[i] * cos(angle) + c[i] * sin(angle);
    }

    for (i = 0; i < DIMENSION; i++){
        to->coords[i] = y[i] + this->distance + this->normal[i];
    }

    to->element = from->element;
    to->pse_num = from->pse_num;

    return 1;
}