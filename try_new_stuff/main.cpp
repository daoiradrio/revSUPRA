#include <iostream>
#include <fstream>
#include <limits>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>



#define TOL 0.000001



void dac(std::vector<double> arr, double search, int l, int r){
        int m;
        while (l <= r){
                m = l + (r - l)/2;
                if (fabs(search - arr[m]) < TOL){
                        std::cout << "Found at " << m << std::endl;
			std::cout << search << " " << arr[m] << std::endl;
			std::cout << std::endl;
                        return;
                }
                else if (search < arr[m]){
                        r = m - 1;
                }
                else{
                        l = m + 1;
                }
        }
        std::cout << "Not found..." << std::endl;
        return;
}



int main(int argc, char **argv){
	int dummy;
	double energy;
	std::string filepath;
	std::string line;
	int i;
	int n_confs = 97;
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

	for (double energy: energies){
		dac(energies, energy, 0, energies.size());
	}

	return 0;
}

