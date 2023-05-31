#include <rotationaxis.hpp>



std::shared_ptr<Atom> RotationAxis::transform_atom(std::shared_ptr<Atom> from){
    int                     i;
    double                  dot = 0.0;
    double                  angle = 0.0;
    std::vector<double>     v1(DIMENSION, 0.0);
    std::vector<double>     v2(DIMENSION, 0.0);
    std::vector<double>     v3(DIMENSION, 0.0);
    std::vector<double>     cross(DIMENSION, 0.0);
    std::shared_ptr<Atom>   to = std::make_shared<Atom>();
    
    if (!this->order){
        std::cout << "Order of rotation axis is zero!" << std::endl;
        return nullptr;
    }

    angle =  (2.0 * M_PI) / (double)this->order;

    Eigen::Vector3d new_coord(from->coords.data());
    Eigen::Vector3d axis_vec1(this->axis_from->coords.data());
    Eigen::Vector3d axis_vec2(this->axis_to->coords.data());
    Eigen::Vector3d axis;
    //Eigen::Vector3d t1;
    //Eigen::Vector3d t2;
    //Eigen::Vector3d t3;
    axis = axis_vec2 - axis_vec1;
    axis.normalize();
    new_coord = new_coord - axis_vec1;
    new_coord = axis.dot(new_coord) * axis
                + cos(angle) * axis.cross(new_coord).cross(axis) // FUNKTIONIERT cross().cross()??
                + sin(angle) * axis.cross(new_coord);
    /*t1 = axis * axis.dot(new_coord);
    t2 = axis.cross(new_coord);
    t2 = t2.cross(axis) * cos(angle);
    t3 = axis.cross(new_coord) * sin(angle);
    new_coord = t1 + t2 + t3;*/
    new_coord = new_coord + axis_vec1;

    /*for (i = 0; i < DIMENSION; i++){
        to->coords[i] = from->coords[i] - this->axis_from->coords[i];
    }

    for (i = 0; i < DIMENSION; i++){
        dot += this->direction[i] * from->coords[i];
    }
    
    for (i = 0; i < DIMENSION; i++){
        v1[i] = this->direction[i] * dot;
    }

    cross[0] = this->direction[1] * from->coords[2] - this->direction[2] * from->coords[1];
    cross[1] = this->direction[2] * from->coords[0] - this->direction[0] * from->coords[2];
    cross[2] = this->direction[0] * from->coords[1] - this->direction[1] * from->coords[0];

    v2[0] = cross[1] * this->direction[2] - cross[2] * this->direction[1];
    v2[1] = cross[2] * this->direction[0] - cross[0] * this->direction[2];
    v2[2] = cross[0] * this->direction[1] - cross[1] * this->direction[0];

    for (i = 0; i < DIMENSION; i++){
        v2[i] *= cos(angle);
    }

    v3[0] = this->direction[1] * from->coords[2] - this->direction[2] * from->coords[1];
    v3[1] = this->direction[2] * from->coords[0] - this->direction[0] * from->coords[2];
    v3[2] = this->direction[0] * from->coords[1] - this->direction[1] * from->coords[0];

    for (i = 0; i < DIMENSION; i++){
        v3[i] *= sin(angle);
    }

    for (i = 0; i < DIMENSION; i++){
        to->coords[i] = v1[i] + v2[i] + v3[i];
    }

    for (i = 0; i < DIMENSION; i++){
        to->coords[i] += this->axis_from->coords[i];
    }*/

    to->element = from->element;
    to->pse_num = from->pse_num;
    for (i = 0; i < DIMENSION; i++){
        to->coords[i] = new_coord(i);
    }

    return to;
}