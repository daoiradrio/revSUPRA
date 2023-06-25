#include <optimizer.hpp>



int Optimizer::uff_optimization(std::string path, std::string xyz_file, int index)
{
    int             fin;
    int             line_index;
    int             dummy;
    double          energy;
    std::string     command;
    std::string     opt_dir;
    std::string     coord_file;
    std::string     control_file;
    std::ifstream   infile;
    std::ofstream   outfile;
    std::string     line;

    if (path.back() != '/'){
        path = path + "/";
    }

    // create working directory for optimization
    if (index != -1){
        opt_dir = path + "struc_opt" + std::to_string(index) + "/";
    }
    else{
        opt_dir = path + "struc_opt/";
    }
    command = "mkdir " + opt_dir;
    fin = system(command.c_str());
    
    if (fin < 0){
        return FAIL_EXIT;
    }

    // move .xyz file of unoptimized target structure to working directory
    //command = "mv " + path + xyz_file + " " + opt_dir;
    //system(command.c_str());

    // convert coordinates to TURBOMOLE format
    coord_file = opt_dir + "coord";
    command = "x2t " + path + xyz_file + " > " + coord_file;
    //command = "x2t " + opt_dir + xyz_file + " > " + coord_file;
    fin = system(command.c_str());

    if (fin < 0){
        return FAIL_EXIT;
    }

    // write control file for UFF optimization
    control_file = opt_dir + "control";
    outfile.open(control_file);
    outfile << "$symmetry c1\n";
    outfile << "$uff\n";
    outfile << "      2500         1          0 ! maxcycle,modus,nqeq\n";
    outfile << "    111111                      ! iterm\n";
    outfile << "  0.10D-07  0.10D-04            ! econv,gconv\n";
    outfile << "      0.00  1.10                ! qtot,dfac\n";
    outfile << "  0.10D+03  0.10D-04       0.30 ! epssteep,epssearch,dqmax\n";
    outfile << "        25      0.10       0.00 ! mxls,dhls,ahls\n";
    outfile << "      1.00      0.00       0.00 ! alpha,beta,gamma\n";
    outfile << "         F         F          F ! transform,lnumhess,lmd\n";
    outfile << "$end\n";
    outfile.close();

    // perform UFF optimization
	command = "cd " + opt_dir + " ; uff > uff.out 2>&1";
    fin = system(command.c_str());

    if (fin < 0){
        return FAIL_EXIT;
    }

    // convert optimized coordinates back to .xyz format and move file to inital directory
    //command = "t2x " + coord_file + " > " + path + xyz_file + " 2>/dev/null";
    //fin = system(command.c_str());
    command = "t2x " + coord_file + " > " + opt_dir + xyz_file + " 2>/dev/null";
    fin = system(command.c_str());

    infile.open(opt_dir + "uffenergy");
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    if (infile.is_open()){
        getline(infile, line);
        std::stringstream linestream(line);
        linestream >> dummy >> energy;
    }
    else{
        return FAIL_EXIT;
    }
    infile.close();

    outfile.open(path + xyz_file);
    infile.open(opt_dir + xyz_file);
    if (infile.is_open() && outfile.is_open()){
        line_index = 0;
        while (getline(infile, line)){
            std::stringstream linestream(line);
            if (line_index == 1){
                outfile << "UFF-Energy = ";
                outfile << energy;
                outfile << "\n";
            }
            else{
                outfile << line;
                outfile << "\n";
            }
            line_index++;
        }
    }
    else{
        return FAIL_EXIT;
    }
    infile.close();
    outfile.close();

    if (fin < 0){
        return FAIL_EXIT;
    }

    // remove working directory
    std::cout << opt_dir << std::endl;
    command = "rm -f " + opt_dir + "*";
    system(command.c_str());
    command = "rm -rf " + opt_dir;
    system(command.c_str());

    if (fin < 0){
        return FAIL_EXIT;
    }

    return SUCCESS_EXIT;
}
