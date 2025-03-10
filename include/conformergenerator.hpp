#ifndef CONFORMERGENERATOR_HPP
#define CONFORMERGENERATOR_HPP

#include <utils.hpp>
#include <structure.hpp>
#include <analyzer.hpp>
#include <rotationaxis.hpp>
#include <symmetry.hpp>
#include <optimizer.hpp>

#include <vector>
#include <memory>
#include <string.h>
#include <cmath>
#include <fstream>
#include <queue>
#include <algorithm>
#include <unistd.h>
#include <Eigen/Dense>



class ConformerGenerator{

    public:
        std::shared_ptr<Structure>          mol;
        Optimizer                           optimizer;
        Analyzer                            analyzer;

        //std::vector<Bond>               torsions;
        //std::vector<Bond>               central_torsions;
        //std::vector<Bond>               methylalike_torsions;
        //std::vector<Bond>               terminal_torsions;

        std::vector<std::shared_ptr<Torsion>>  torsions;
        std::vector<std::shared_ptr<Torsion>>  central_torsions;
        std::vector<std::shared_ptr<Torsion>>  methylalike_torsions;
        std::vector<std::shared_ptr<Torsion>>  terminal_torsions;

        std::vector<std::vector<int>>   torsion_atoms;

        std::vector<int>                angle_increments;
        std::vector<int>                angles;
        std::vector<std::vector<int>>   rot_angles;

        std::vector<Eigen::Vector3d>    input_coords;
        Eigen::MatrixX3d                input_coords_mat;

        std::string                     curr_work_dir;
        std::string                     workdir_name = "opt_dir";
        std::string                     struc_filename = "conformer";
        std::string                     output_foldername = "SUPRA_Output";
        std::string                     opt_struc_filename = "opt_struc.xyz";

        ConformerGenerator(Structure input_mol);
        ~ConformerGenerator();

        void                generate_conformers();
        void                get_torsions();
        void                find_cycles();
        void                cycle_detection(int current, int last, char status[], int ancestors[]);
        void                find_peptidebonds();
        void                selection_menu();
        void                generation_setup();
        std::vector<int>    torsion_atom_counter(int start, int last, std::vector<int> status, std::vector<int> container);
	    std::vector<int>    get_torsion_group(int start, int last, std::vector<int> status, std::vector<int> container = {});
        void                check_rot_sym(int angle_increment);
        int                 combinations(std::vector<Eigen::Vector3d> new_coords, int index, int counter);
        int                 combinations(Eigen::MatrixX3d new_coords, int index, int counter);
        int                 old_combinations(Eigen::MatrixX3d new_coords, int index, int counter);
        bool                clashes(std::vector<Eigen::Vector3d> coords);
        bool                clashes(Eigen::MatrixX3d coords);
        bool                distant_atoms(int atom1, int atom2);
        void                write_xyz(std::vector<Eigen::Vector3d> coords, std::string destination);
        void                write_xyz(Eigen::MatrixX3d coords, std::string destination);
    
    private:
        /*std::shared_ptr<Structure> mol;

        std::vector<Bond> torsions;
        std::vector<Bond> central_torsions;
        std::vector<Bond> methylalike_torsions;
        std::vector<Bond> terminal_torsions;
        std::vector<std::vector<int>> torsion_atoms;

        std::vector<double> angle_increments; // IN RAD STATT DEGREE??
        std::vector<double> angles;

        std::vector<Eigen::Vector3d> input_coords;

        std::string workdir_name;
        std::string struc_filename;
        std::string opt_struc_filename;
        //std::string coord_file;
        //std::string control_file;

        void get_torsions();

        void find_cycles();
        void cycle_detection(int current, int last, char status[], int ancestors[]);

        void find_peptidebonds();

        void selection_menu();
        void generation_setup();
        std::vector<int> torsion_atom_counter(int start, int last, int* status, std::vector<int> container);

        int combinations(std::vector<Eigen::Vector3d> new_coords, int index, int counter);

        bool clashes(std::vector<Eigen::Vector3d> coords);
        bool distant_atoms(int atom1, int atom2);

        void write_xyz(std::vector<Eigen::Vector3d> coords, std::string destination);

        void uff_optimization(std::string workdir);*/
};


#endif
