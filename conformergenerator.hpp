#include "structure.cpp"



class ConformerGenerator{

    public:
        ConformerGenerator();
        ~ConformerGenerator();

        std::vector<bond> central_torsions;
        std::vector<bond> methyl_torsions;
        std::vector<bond> terminal_torsions;

        void generate_conformers(Structure* molecule);
    
    private:
        std::vector<bond> torsions;

        void get_torsions_based_on_bond_matrix();
        void get_torsions(std::vector<bond>* bonds);
};