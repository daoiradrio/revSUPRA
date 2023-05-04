#include "structure.cpp"



class ConformerGenerator{

    public:
        //ConformerGenerator(std::shared_ptr<Structure> input_mol);
        //ConformerGenerator(Structure* input_mol);
        ConformerGenerator(Structure input_mol);
        ~ConformerGenerator();

        void generate_conformers();

        std::vector<bond> ring_bonds;

        std::vector<bond> central_torsions;
        std::vector<bond> methylalike_torsions;
        std::vector<bond> terminal_torsions;

        std::shared_ptr<Structure> mol;
        //Structure* mol;
    
    private:
        //std::shared_ptr<Structure> mol;
        std::vector<bond> torsions;

        void get_torsions();
        void find_cycles();
        void cycle_detection(int current, int last, char status[], int ancestors[]);
};