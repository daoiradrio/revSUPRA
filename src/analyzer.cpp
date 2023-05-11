#include <analyzer.hpp>


Analyzer::Analyzer(){};


Analyzer::~Analyzer(){};


void Analyzer::read_xyz(std::string filepath){
    std::string element;
    std::string line;
    std::ifstream file(filepath);
    float xcoord, ycoord, zcoord;

    if (file.is_open()){
        getline(file, line);
        getline(file, line);
        while (getline(file, line)){
            std::stringstream linestream(line);
            linestream >> element >> xcoord >> ycoord >> zcoord;
            this->elements.push_back(element);
            this->coords.push_back({xcoord, ycoord, zcoord});
        }
    }
    else{
        std::cout << "FAILED OPENING .xyz FILE!" << std::endl;
    }
    file.close();

    return;
}


void Analyzer::kabsch(){
    // ... 
    return;
}


float Analyzer::rmsd(){
    // ...
    return 1.0;
}
