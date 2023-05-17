#include <iostream>
#include <fstream>
#include <limits>
#include <sstream>
#include <vector>
#include <algorithm>



int main(int argc, char **argv){
	int dummy;
	double energy;
	std::string filepath;
	std::string line;
	int i;
	int n_confs = 240;
	std::vector<double> energies;
	
	for (i = 0; i < n_confs; i++){
		filepath = "opt_dir" + std::to_string(i) + "/uffenergy";
		std::ifstream file(filepath);
		file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		if (file.is_open()){
			getline(file, line);
			std::stringstream linestream(line);
			linestream >> dummy >> energy;
			energies.push_back(energy);
		}
	}

	for (i = 0; i < energies.size(); i++){
		std::cout << energies[i] << std::endl;
	}
	std::cout << std::endl;

	std::sort(energies.begin(), energies.end());

	for (i = 0; i < energies.size(); i++){
                std::cout << energies[i] << std::endl;
        }
	std::cout << std::endl;

	return 0;
}

