#include <rotationaxis.hpp>



RotationAxis::RotationAxis(std::shared_ptr<Atom> from_atom, std::shared_ptr<Atom> to_atom){
    this->from << from_atom->coords[0], from_atom->coords[1], from_atom->coords[2];
    this->to << to_atom->coords[0], to_atom->coords[1], to_atom->coords[2];
    this->axis = this->to - this->from;
    this->axis.normalize();
}



RotationAxis::RotationAxis(std::vector<double> from_coords, std::vector<double> to_coords){
    this->from << from_coords[0], from_coords[1], from_coords[2];
    this->to << to_coords[0], to_coords[1], to_coords[2];
    this->axis = this->to - this->from;
    this->axis.normalize();
}



RotationAxis::RotationAxis(Eigen::Vector3d from_coords, Eigen::Vector3d to_coords){
    this->from = from_coords;
    this->to = to_coords;
    this->axis = this->to - this->from;
    this->axis.normalize();
}



Eigen::Vector3d RotationAxis::rotate_atom(Eigen::Vector3d coords, double deg){
    Eigen::Vector3d     new_coords;
    double              angle = (2.0 * M_PI) * (deg / 360.0);
    
    new_coords = coords - this->from;
    new_coords = axis.dot(new_coords) * axis
                + cos(angle) * this->axis.cross(new_coords).cross(this->axis)
                + sin(angle) * this->axis.cross(new_coords);
    new_coords = new_coords + this->from;

    return new_coords;
}



std::vector<double> RotationAxis::rotate_atom(std::vector<double> coords, double deg){
    std::vector<double>     new_coords(DIMENSION, 0.0);
    Eigen::Vector3d         parsed_coords(coords.data());

    parsed_coords = this->rotate_atom(parsed_coords, deg);
    new_coords[0] = parsed_coords(0);
    new_coords[1] = parsed_coords(1);
    new_coords[2] = parsed_coords(2);

    return new_coords;
}



std::shared_ptr<Atom> RotationAxis::rotate_atom(std::shared_ptr<Atom> atom, double deg){
    std::shared_ptr<Atom>   new_atom = std::make_shared<Atom>();
    Eigen::Vector3d         new_coords(atom->coords.data());
    
    new_coords = this->rotate_atom(new_coords, deg);
    new_atom->coords[0] = new_coords(0);
    new_atom->coords[1] = new_coords(1);
    new_atom->coords[2] = new_coords(2);
    new_atom->element = atom->element;
    new_atom->pse_num = atom->pse_num;

    return new_atom;
}
