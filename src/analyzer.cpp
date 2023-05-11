#include <analyzer.hpp>


Analyzer::Analyzer(){};


Analyzer::~Analyzer(){};


void Analyzer::read_xyz(std::string filepath){
    std::string element;
    std::string line;
    std::ifstream file(filepath);
    double xcoord, ycoord, zcoord;
    Eigen::Vector3d new_coords;

    if (file.is_open()){
        getline(file, line);
        getline(file, line);
        while (getline(file, line)){
            std::stringstream linestream(line);
            linestream >> element >> xcoord >> ycoord >> zcoord;
            new_coords = {xcoord, ycoord, zcoord};
            this->elements.push_back(element);
            //this->coords.push_back(new_coords);
        }
    }
    else{
        std::cout << "FAILED OPENING .xyz FILE!" << std::endl;
    }
    file.close();

    return;
}


void Analyzer::remove_doubles(){
    return;
}


void Analyzer::kabsch(Eigen::Matrix3d coords1, Eigen::Matrix3d coords2){
    int i;
    Eigen::Vector3d center1 = {0.0, 0.0, 0.0};
    Eigen::Vector3d center2 = {0.0, 0.0, 0.0};

    /*
    for (i = 0; i < coords1.size(); i++){
        center1 = center1 + coords1[i];
        center2 = center2 + coords2[i];
    }
    */

    center1 = (1.0/(double)coords1.size()) * center1;
    center2 = (1.0/(double)coords2.size()) * center2;

    /*
    for (i = 0; i < coords1.size(); i++){
        coords1[i] = coords1[i] - center1;
        coords2[i] = coords2[i] - center2;
    }
    */

    return;
}


float Analyzer::rmsd(){
    // ...
    return 1.0;
}
