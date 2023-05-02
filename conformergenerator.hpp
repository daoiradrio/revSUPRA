#include "structure.cpp"



class ConformerGenerator{

    public:
        ConformerGenerator();
        ~ConformerGenerator();

        void generate_conformers(Structure* molecule);

        std::vector<bond*> ring_bonds;

        std::vector<bond*> central_torsions;
        std::vector<bond*> methylalike_torsions;
        std::vector<bond*> terminal_torsions;
    
    private:
        std::vector<bond> torsions;

        void get_torsions(std::vector<bond> bonds);
        //void find_cycles(std::vector<atom*> atoms);
        //void cycle_detection(std::vector<atom*> atoms, atom* current, atom* last, char status[], int ancestors[]);
};