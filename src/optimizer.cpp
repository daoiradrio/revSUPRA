#include <optimizer.hpp>



void Optimizer::uff_optimization(std::string curr_dir, std::string xyz_file)
{
    int             fin;
    std::string     command;
    std::string     workdir;
    std::string     xyz_file_path;
    std::string     coord_file;
    std::string     control_file;

    if (curr_dir.back() != '/'){
        curr_dir = curr_dir + "/";
    }

    workdir = curr_dir + this->workdir_name;
    command = "mkdir " + workdir;
    system(command.c_str());

    command = "mv " + curr_dir + xyz_file + " " + workdir;
    system(command.c_str());

    coord_file = workdir + "coord";
    command = "x2t " + xyz_file_path + " > " + coord_file;
    system(command.c_str());

    // write control file for UFF optimization
    control_file = workdir + "control";
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
	command = "cd " + workdir + " ; uff > uff.out 2>&1 &";
    fin = system(command.c_str());

    command = "t2x " + coord_file + " > " + curr_dir + xyz_file + " 2>/dev/null";
    system(command.c_str());

    command = "rm -r " + workdir;
    system(command.c_str());

    return;
}
