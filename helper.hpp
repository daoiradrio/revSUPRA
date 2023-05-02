#include <iostream>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <cmath>
#include <sstream>
#include <string.h>



std::unordered_map<std::string, float> valence_radii_single({
    {"C", 0.75},
    {"N", 0.71},
    {"O", 0.66},
    {"H", 0.32},
    {"B", 0.85},
    {"F", 0.64},
    {"Cl", 0.99},
    {"Br", 1.14},
    {"I", 1.33}
});


std::unordered_map<std::string, float> valence_radii_double({
    {"C", 0.67},
    {"N", 0.60},
    {"O", 0.57},
    {"H", NULL},
    {"B", 0.78},
    {"F", NULL},
    {"Cl", NULL},
    {"Br", NULL},
    {"I", NULL}
});


std::unordered_map<std::string, float> valence_radii_triple({
    {"C", 0.60},
    {"N", 0.54},
    {"O", 0.53},
    {"H", NULL},
    {"B", 0.73},
    {"F", NULL},
    {"Cl", NULL},
    {"Br", NULL},
    {"I", NULL}
});


std::unordered_map<std::string, int> max_valences({
    {"C", 4},
    {"N", 3},
    {"O", 2},
    {"B", 3},
    {"H", 1},
    {"F", 1},
    {"Cl", 1},
    {"Br", 1},
    {"I", 1},
});


std::string get_element(std::string label)
{   
    std::string element;
    
    if (isalpha(label[1])){
        element = std::string() + label[0] + label[1];
    }
    else{
        element = label[0];
    }

    return element;
}


int hash_bond_matrix(int row_index, int column_index, int n_atoms){
    return (2*row_index*n_atoms - row_index - pow(row_index, 2))/2 + column_index - row_index - 1;
}


bool is_terminal_atom(std::string element){
    bool terminal = false;
    if (element == "H"){
        terminal = true;
    }
    else if (element == "F"){
        terminal = true;
    }
    else if (element == "Cl"){
        terminal = true;
    }
    else if (element == "Br"){
        terminal = true;
    }
    else if (element == "I"){
        terminal = true;
    }
    return terminal;
}


typedef struct atom{
    std::string element;
    int index;
    std::string label;
    double coords[3];
    std::vector<atom> bond_partners;
    bool core_of_terminal_group = false;
}atom;


typedef struct bond{
    atom atom1;
    atom atom2;
    int bond_order;
}bond;