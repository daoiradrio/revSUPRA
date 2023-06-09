#include <optimizer.hpp>



int Optimizer::uff_optimization(std::string path, std::string xyz_file, int index)
{
    int             fin;
    std::string     command;
    std::string     opt_dir;
    std::string     coord_file;
    std::string     control_file;
    std::ifstream   output_file;

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
	std::ofstream file;
    file.open(control_file);
    file << "$symmetry c1\n";
    file << "$uff\n";
    file << "      2500         1          0 ! maxcycle,modus,nqeq\n";
    file << "    111111                      ! iterm\n";
    file << "  0.10D-07  0.10D-04            ! econv,gconv\n";
    file << "      0.00  1.10                ! qtot,dfac\n";
    file << "  0.10D+03  0.10D-04       0.30 ! epssteep,epssearch,dqmax\n";
    file << "        25      0.10       0.00 ! mxls,dhls,ahls\n";
    file << "      1.00      0.00       0.00 ! alpha,beta,gamma\n";
    file << "         F         F          F ! transform,lnumhess,lmd\n";
    file << "$end\n";
    file.close();

    // perform UFF optimization
	command = "cd " + opt_dir + " ; uff > uff.out 2>&1";
    fin = system(command.c_str());

    /*output_file.open(opt_dir + "uff.out");
    while (!output_file.is_open()){
        continue;
    }
    output_file.close();*/

    if (fin < 0){
        return FAIL_EXIT;
    }

    // convert optimized coordinates back to .xyz format and move file to inital directory
    command = "t2x " + coord_file + " > " + path + xyz_file + " 2>/dev/null";
    //command = "t2x " + coord_file + " > " + opt_dir + "opt_struc.xyz" + " 2>/dev/null";
    fin = system(command.c_str());

    if (fin < 0){
        return FAIL_EXIT;
    }

    // remove working directory
    command = "rm -f " + opt_dir + "*";
    system(command.c_str());
    command = "rm -rf " + opt_dir;
    system(command.c_str());

    if (fin < 0){
        return FAIL_EXIT;
    }

    return SUCCESS_EXIT;
}
