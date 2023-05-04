#include "structure.cpp"



class ConformerGenerator{

    public:
        ConformerGenerator(std::unique_ptr<Structure>* mol);
        ~ConformerGenerator();

        void generate_conformers();

        std::vector<bond> ring_bonds;

        std::vector<bond> central_torsions;
        std::vector<bond> methylalike_torsions;
        std::vector<bond> terminal_torsions;
    
    private:
        std::unique_ptr<Structure> mol;
        std::vector<bond> torsions;

        void get_torsions();
        //void find_cycles(std::vector<atom*> atoms);
        //void cycle_detection(std::vector<atom*> atoms, atom* current, atom* last, char status[], int ancestors[]);
};