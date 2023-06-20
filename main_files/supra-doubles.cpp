#include <analyzer.hpp>

#include <iostream>



void print_usage()
{   
    std::cout                                                                             << std::endl;
    std::cout << "USAGE:"                                                                 << std::endl;
    std::cout << "supra-doubles -path /path/to/directory"                                 << std::endl;
    std::cout << "supra-doubles -path /path/to/directory -ignore methyl"                  << std::endl;
    std::cout << "supra-doubles -path1 /path/to/file -path2 /path/to/file"                << std::endl;
    std::cout << "supra-doubles -path1 /path/to/file -path2 /path/to/file -ignore methyl" << std::endl;
    std::cout                                                                             << std::endl;

    return;
}



int main(int argc, char** argv)
{   
    // supra-doubles -path PATH
    // supra-doubles -path1 PATH -path2 PATH
    // supra-doubles -path PATH -ignore methyl
    // supra-doubles -path1 PATH -path2 PATH -ignore methyl
    // supra-doubles -path PATH -rmsd X
    // supra-doubles -path1 PATH -path2 PATH -rmsd X
    // supra-doubles -path PATH -ignore methyl -rmsd X
    // supra-doubles -path1 PATH -path2 PATH -ignore methyl -rmsd X

    Analyzer analyzer;

    int         i;
    std::string mode; // "structure_pair" or "structure_set"
    std::string path1;
    std::string path2;
    bool        ignore_methyl = false;
    double      rmsd;
    double      rmsd_threshold = 0.1;

    if (argc < 3){
        print_usage();
        return 0;
    }

    i = 1;
    while (i < argc){
        if (std::string(argv[i]) == "-path"){
            mode = "structure_set";
            i++;
            if (i < argc){
                path1 = std::string(argv[i]);
                i++;
            }
            else{
                print_usage();
                return 0;
            }
        }
        else if (std::string(argv[i]) == "-path1"){
            mode = "structure_pair";
            i++;
            if (i < argc){
                path1 = std::string(argv[i]);
                i++;
            }
            else{
                print_usage();
                return 0;
            }
        }
        else if (std::string(argv[i]) == "-path2"){
            mode = "structure_pair";
            i++;
            if (i < argc){
                path2 = std::string(argv[i]);
                i++;
            }
            else{
                print_usage();
                return 0;
            }
        }
        else if (std::string(argv[i]) == "-ignore"){
            i++;
            if (i < argc && std::string(argv[i]) == "methyl"){
                ignore_methyl = true;
                i++;
            }
            else{
                print_usage();
                return 0;
            }
        }
        else if (std::string(argv[i]) == "-rmsd"){
            i++;
            if (i < argc){
                rmsd_threshold = std::stod(argv[i]);
                i++;
            }
            else{
                print_usage();
                return 0;
            }
        }
        else{
            print_usage();
            return 0;
        }
    }

    if (mode == "structure_pair"){
        if (ignore_methyl){
            if (analyzer.doubles(path1, path2, true, rmsd_threshold)){
                std::cout << "Doubles" << std::endl;
            }
        }
        else{
            if (analyzer.doubles(path1, path2, rmsd_threshold)){
                std::cout << "No doubles" << std::endl;
            }
        }
    }
    else if (mode == "structure_set"){
        analyzer.remove_doubles(path1, "conformer", ignore_methyl, rmsd_threshold);
    }
    else{
        print_usage();
    }
    
    return 0;
}