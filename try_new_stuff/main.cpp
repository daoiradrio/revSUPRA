#include <iostream>
#include <fstream>
#include <limits>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>



#define TOL 0.00005
#define EQ_DBL_DELTA 1E-16



struct cmp{
    cmp(std::vector<double> vec) : val(vec[0]){}
    inline bool operator()(std::vector<double>& x){
        return (fabs(x[0]-val) < EQ_DBL_DELTA);
    }
private:
    double val;
};



int dac(std::vector<double> arr, double search, int l, int r){
    int m;
    while (l <= r){
        m = l + (r - l)/2;
        if (fabs(search - arr[m]) < TOL){
            return 1;
        }
        else if (search < arr[m]){
            r = m - 1;
        }
        else{
            l = m + 1;
        }
    }
    return 0;
}



int binary_search(std::vector<std::vector<double>> container, double val, int l, int r){
    int m;
    while (l <= r){
        m = l + (r - l)/2;
        if (fabs(val - container[m][0]) < TOL){
            //return container[m][1];
            return 1;
        }
        else if (val < container[m][0]){
            r = m - 1;
        }
        else{
            l = m + 1;
        }
    }
    //return -1;
    return 0;
}



bool sort_func(std::vector<double> a, std::vector<double> b){
    return (a[0] < b[0]);
}



int main(int argc, char **argv){
	int i, j;
	int dummy;
	int n_confs = 97;
    double energy;
	std::string filepath;
	std::string line;
	std::vector<double> energies;
    std::vector<double> copy_energies;
    std::vector<std::vector<double>> container;
	
	for (i = 0; i < n_confs; i++){
		filepath = "opt_dir" + std::to_string(i) + "/uffenergy";
		std::ifstream file(filepath);
		file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		if (file.is_open()){
			getline(file, line);
			std::stringstream linestream(line);
			linestream >> dummy >> energy;
			energies.push_back(energy);
            container.push_back({energy, (double)i});
		}
	}

	std::sort(energies.begin(), energies.end());
    std::sort(container.begin(), container.end(), sort_func);

    int counter = 0;
    std::vector<double>::iterator it1;
    for (double energy: energies){
        copy_energies = energies;
        //it = std::find_if(copy_energies.begin(), copy_energies.end(), dbl_cmp(energy));
        it1 = std::find(copy_energies.begin(), copy_energies.end(), energy);
        copy_energies.erase(it1);
        if (dac(copy_energies, energy, 0, copy_energies.size())){
            counter++;
        }
    }
    std::cout << "Von " << n_confs << " sind " << counter << " Doubles." << std::endl;
    counter = 0;
    std::vector<std::vector<double>> copy_container;
    std::vector<std::vector<double>>::iterator it2;
    std::vector<double> item;
    for (i = 0; i < n_confs; i++){
        copy_container = container;
        item = copy_container[i];
        energy = item[0];
        it2 = std::find_if(copy_container.begin(), copy_container.end(), cmp(item));
        copy_container.erase(it2);
        if (binary_search(copy_container, energy, 0, copy_container.size()-1)){
            counter++;
        }
    }
    std::cout << "Von " << n_confs << " sind " << counter << " Doubles." << std::endl;

	return 0;
}

