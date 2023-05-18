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
    int i;
	int n_confs = 97;
    double energy;
	std::vector<double> energies;
    std::string path;
    std::string file;
    std::string line;
    std::string bash_script = "extract_energy.sh";
    std::ofstream bash_script_stream;
    std::ifstream energy_file_stream;
    std::string command;

    bash_script_stream.open(bash_script);
    bash_script_stream << "#!/bin/bash\n";
    bash_script_stream << "gawk '/TOTAL ENERGY/ { val=$4; printf(\"%18.12f\", val); }' $1 >> $2";
    bash_script_stream.close();
    command = "chmod +x " + bash_script;
    system(command.c_str());
	
	for (i = 0; i < n_confs; i++){
        path = "opt_dir" + std::to_string(i) + "/gfnff_singlepoint/";
        file = path + "gfnff.out";
        command = "./" + bash_script + " " + file + " " + path + ".EOUT";
        system(command.c_str());
        file = path + ".EOUT";
        command = "cat " + file;
        system(command.c_str());
        energy_file_stream.open(file);
        if (energy_file_stream.is_open()){
            getline(energy_file_stream, line);
            std::stringstream linestream(line);
            linestream >> energy;
            energies.push_back(energy);
        }
	}

	//for (i = 0; i < energies.size(); i++){
	//	std::cout << energies[i] << std::endl;
	//}
	std::cout << std::endl;

	//std::sort(energies.begin(), energies.end());

	//for (i = 0; i < energies.size(); i++){
    //            std::cout << energies[i] << std::endl;
    //    }
	//std::cout << std::endl;

	//for (double energy: energies){
	//	dac(energies, energy, 0, energies.size());
	//}

	return 0;
}

