#include <symmetry.hpp>



Symmetry::Symmetry(){}



Symmetry::~Symmetry(){}



int Symmetry::establish_pairs(std::shared_ptr<SYMMETRY_ELEMENT> elem)
{
    int                     i, j, k;
    int                     best_j;
    double                  distance, best_distance;
    std::vector<int>        atom_used(this->AtomsCount, 1);
    std::shared_ptr<ATOM>   symmetric;

    for (i = 0; i < this->AtomsCount; i++){
        if (elem->transform[i] >= this->AtomsCount){
            continue;
        }
        elem->transform_atom(elem, this->atoms[i], symmetric); // HIER TESTEN OB BEIDE VERSIONEN GLEICH LAUFEN
    }
}



void Symmetry::find_geometric_center()
{
    int                 i, j;
    double              r;
    std::vector<double> coord_sum(DIMENSION, 0.0);

    for (i = 0; i < this->AtomsCount; i++){
        for (j = 0; j < DIMENSION; j++){
            coord_sum[j] += this->atoms[i]->coords[j];
        }
    }

    for (i = 0; i < DIMENSION; i++){
        this->geom_center[i] = coord_sum[i] / this->AtomsCount;
    }

    for (i = 0; i < this->AtomsCount; i++){
        r = 0.0;
        for (j = 0; j < DIMENSION; j++){
            r += pow(this->atoms[i]->coords[j] - this->geom_center[j], 2);
        }
        this->dist_geom_center[i] = r;
    }

    return;
}



void Symmetry::check_C2_axis()
{
    int                                 i, j, k;
    int                                 sum;
    double                              distance_i, distance_j;
    std::vector<int>                    flags(this->AtomsCount, 0);
    std::shared_ptr<SYMMETRY_ELEMENT>   axis;

    for (i = 1; i < this->AtomsCount; i++){
        for (j = 0; j < i; j++){
            if (this->atoms[i]->pse_num != this->atoms[j]->pse_num){
                continue;
            }
            if (fabs(this->dist_geom_center[i] - this->dist_geom_center[j]) > TolerancePrimary){
                continue;
            }
            distance_i = 0.0;
            distance_j = 0.0;
            for (k = 0; k < DIMENSION; k++){
                distance_i += pow(this->atoms[i]->coords[k] - this->support[k], 2);
                distance_j += pow(this->atoms[j]->coords[k] - this->support[k], 2);
            }
            distance_i = sqrt(distance_i);
            distance_j = sqrt(distance_j);
            if ((distance_i - distance_j) > TolerancePrimary){
                continue;
            }
            if (this->init_C2(i, j, support)){
                flags[i] = 1;
                flags[j] = 1;
            }
        }
    }

    return;
}



int Symmetry::init_C2(int atom1, int atom2, std::vector<double> support)
{
    std::shared_ptr<SYMMETRY_ELEMENT>   axis;
    int                                 i;
    double                              r, r1, r2;
    std::vector<double>                 middle(DIMENSION, 0.0);

    std::fill(axis->transform.begin(), axis->transform.end(), this->AtomsCount+1);

    r1 = 0.0;
    r2 = 0.0;
    for (i = 0; i < DIMENSION; i++){
        r1 += pow(this->atoms[atom1]->coords[i] - support[i], 2);
        r2 += pow(this->atoms[atom2]->coords[i] - support[i], 2);
    }
    r1 = sqrt(r1);
    r2 = sqrt(r2);
    if (fabs(r1 - r2) > TolerancePrimary){
        return 0;
    }

    r = 0.0;

    for (i = 0; i < DIMENSION; i++){
        r += pow(this->geom_center[i], 2);
    }
    r = sqrt(r);
    if (r > 0.0){
        for (i = 0; i < DIMENSION; i++){
            axis->normal[i] = this->geom_center[i] / r;
        }
    }
    else{
        axis->normal[0] = 1.0;
        axis->normal[1] = 0.0;
        axis->normal[2] = 0.0;
    }

    axis->distance = r;
    r = 0.0;

    for (i = 0; i < DIMENSION; i++){
        middle[i] =  (this->atoms[atom1]->coords[i] + this->atoms[atom2]->coords[i]) / 2 - support[i];
        r         += middle[i] * middle[i];
    }
    r = sqrt(r);

    if (r <= TolerancePrimary){
        for (i = 0; i < DIMENSION; i++){
            middle[i] = this->atoms[atom1]->coords[i] - this->atoms[atom2]->coords[i];
        }
        if (fabs(middle[2]) + fabs(middle[1]) > ToleranceSame){
            axis->direction[0] =  0.0;
            axis->direction[1] =  middle[2];
            axis->direction[2] = -middle[1];
        }
        else{
            axis->direction[0] = -middle[2];
            axis->direction[1] =  0.0;;
            axis->direction[2] = -middle[0];
        }
        r = 0.0;
        for (i = 0; i < DIMENSION; i++){
            r += axis->direction[i] * axis->direction[i];
        }
        r = sqrt(r);
        for (i = 0; i < DIMENSION; i++){
            axis->direction[i] /= r;
        }
    }
    else{
        for (i = 0; i < DIMENSION; i++){
            axis->direction[i] = middle[i] / r;
        }
    }

    this->establish_pairs(axis);
    
    return 1;
}



bool Symmetry::detect_rot_sym(std::shared_ptr<Structure> mol, std::vector<int> torsion_atoms, int order) // ANGLE STATT ORDER
{   
    this->support_atom = torsion_atoms[0];
    this->AtomsCount = torsion_atoms.size();

    for (int atom: torsion_atoms){
        this->atoms.push_back(std::make_shared<ATOM>(*(mol->atoms[atom])));
    }

    this->find_geometric_center();

    if (order == 2){
        this->check_C2_axis();
    }
    //else{
    //    this->Cn_axis();
    //}

    return false;
}