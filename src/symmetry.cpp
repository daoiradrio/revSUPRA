#include <symmetry.hpp>



Symmetry::Symmetry(){}



Symmetry::~Symmetry(){}



int Symmetry::establish_pairs(std::shared_ptr<SymmetryElement> elem)
{
    int                     i, j, k;
    int                     best_j;
    double                  distance, best_distance;
    std::vector<int>        atom_used(this->AtomsCount, 1);
    std::shared_ptr<Atom>   symmetric;

    for (i = 0; i < this->AtomsCount; i++){
        if (elem->transform[i] >= this->AtomsCount){
            continue;
        }
        elem->transform_atom(this->atoms[i], symmetric);
        best_j = i;
        best_distance = 2 * TolerancePrimary;
        for (j = 0; j < this->AtomsCount; j++){
            if (this->atoms[i]->pse_num != symmetric->pse_num || atom_used[j]){
                continue;
            }
            distance = 0.0;
            for (k = 0; k < DIMENSION; k++){
                distance = distance + pow(symmetric->coords[k] - this->atoms[k]->coords[k], 2);
            }
            distance = sqrt(distance);
            if (distance < best_distance){
                best_j = j;
                best_distance = distance;
            }
        }
        if (best_distance > TolerancePrimary){
            return 0;
        }
        elem->transform[i] = j;
        atom_used[i] = 1;
    }
    
    return 1;
}



int Symmetry::check_transform_order(std::shared_ptr<SymmetryElement> elem){
    int i, j, k;

    for (i = 0; i < this->AtomsCount; i++){
        if (elem->transform[i] == i){
            continue;
        }
        for (j = elem->order-1, k = elem->transform[i]; j > 0; j--, k = elem->transform[k]){
            if (k == i){
                return 0;
            }
        }
        if (k != i){
            return 0;
        }
    }

    return 1;
}



int Symmetry::optimize_transform_params(std::shared_ptr<SymmetryElement> elem){
    int i;
    int cycle = 0;
    int hits = 0;
    int finish = 0;
    double f = 0;
    double fold, snorm;

    f = 0.0;
    do{
        fold = f;
        f    = this->eval_opt_target_func(elem, std::make_shared<int>(finish));
        if (finish){
            break;
        }
        if (cycle > 0){
            if (fabs(f - fold) > OptChangeThreshold){
                hits = 0;
            }
            else{
                hits++;
            }
            if (hits >= OptChangeHits){
                break;
            }
        }
        
    }while(snorm > MinOptStep && cycle++ < MaxOptCycles);

    // DUMMY
    return 0;
}



double Symmetry::eval_opt_target_func(std::shared_ptr<SymmetryElement> elem, std::shared_ptr<int> finish){
    int                     i, j, k;
    double                  target, r, maxr;
    std::shared_ptr<Atom>   symmetric;

    r = 0.0;
    for (i = 0; i < DIMENSION; i++){
        r = r + pow(elem->normal[i], 2);
    }
    r = sqrt(r);
    for (i = 0; i < DIMENSION; i++){
        elem->normal[i] = elem->normal[i] / r;
    }
    if (elem->distance < 0){
        elem->distance = -elem->distance;
        for (i = 0; i < DIMENSION; i++){
            elem->normal[i] = -elem->normal[i];
        }
    }

    target = 0.0;
    maxr   = 0.0;
    for (i = 0; i < this->AtomsCount; i++){
        elem->transform_atom(this->atoms[i], symmetric);
        j = elem->transform[i];
        r = 0.0;
        for(j = 0; j < DIMENSION; j++){
            r = r + pow(this->atoms[j]->coords[j] - symmetric->coords[j], 2);
        }
        if (r > maxr){
            maxr = r;
        }
        target = target + r;
        *finish = 0;
        if (sqrt(maxr) < ToleranceFinal){
            *finish = 1;
        }
    }

    return target;
}



void Symmetry::get_params(std::shared_ptr<SymmetryElement> elem, std::vector<double> values){
    return;
}



void Symmetry::set_params(std::shared_ptr<SymmetryElement> elem, std::vector<double> values){
    return;
}



void Symmetry::find_geometric_center()
{
    int                 i, j;
    double              r;
    std::vector<double> coord_sum(DIMENSION, 0.0);

    this->geom_center.resize(DIMENSION);
    this->dist_geom_center.resize(this->AtomsCount);

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
    //std::shared_ptr<SYMMETRY_ELEMENT>   axis;
    //RotationAxis                        axis;
    //std::shared_ptr<RotationAxis>       axis_ptr = std::make_shared<RotationAxis>(axis);

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
            if (this->init_C2(i, j)){
                flags[i] = 1;
                flags[j] = 1;
            }
        }
    }

    return;
}



// MUSS AXIS TATSÄCHLICH ALS POINTER VERWENDET WERDEN??
int Symmetry::init_C2(int atom1, int atom2)
{   
    //SYMMETRY_ELEMENT                    axis;
    //std::shared_ptr<SYMMETRY_ELEMENT>   axis_ptr = std::make_shared<SYMMETRY_ELEMENT>(axis);
    RotationAxis                        axis;
    std::shared_ptr<RotationAxis>       axis_ptr = std::make_shared<RotationAxis>(axis);
    int                                 i;
    double                              r, r1, r2;
    std::vector<double>                 middle(DIMENSION, 0.0);

    axis_ptr->order = 2;
    axis_ptr->transform.resize(this->AtomsCount);
    std::fill(axis_ptr->transform.begin(), axis_ptr->transform.end(), this->AtomsCount+1);
    //axis_ptr->transform_atom = &rotate_atom;

    r1 = 0.0;
    r2 = 0.0;
    for (i = 0; i < DIMENSION; i++){
        r1 += pow(this->atoms[atom1]->coords[i] - this->support[i], 2);
        r2 += pow(this->atoms[atom2]->coords[i] - this->support[i], 2);
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
            axis_ptr->normal[i] = this->geom_center[i] / r;
        }
    }
    else{
        axis_ptr->normal[0] = 1.0;
        axis_ptr->normal[1] = 0.0;
        axis_ptr->normal[2] = 0.0;
    }

    axis_ptr->distance = r;
    r = 0.0;

    for (i = 0; i < DIMENSION; i++){
        middle[i] =  (this->atoms[atom1]->coords[i] + this->atoms[atom2]->coords[i]) / 2 - this->support[i];
        r         += middle[i] * middle[i];
    }
    r = sqrt(r);

    if (r <= TolerancePrimary){
        for (i = 0; i < DIMENSION; i++){
            middle[i] = this->atoms[atom1]->coords[i] - this->atoms[atom2]->coords[i];
        }
        if (fabs(middle[2]) + fabs(middle[1]) > ToleranceSame){
            axis_ptr->direction[0] =  0.0;
            axis_ptr->direction[1] =  middle[2];
            axis_ptr->direction[2] = -middle[1];
        }
        else{
            axis_ptr->direction[0] = -middle[2];
            axis_ptr->direction[1] =  0.0;;
            axis_ptr->direction[2] = -middle[0];
        }
        r = 0.0;
        for (i = 0; i < DIMENSION; i++){
            r += axis_ptr->direction[i] * axis_ptr->direction[i];
        }
        r = sqrt(r);
        for (i = 0; i < DIMENSION; i++){
            axis_ptr->direction[i] /= r;
        }
    }
    else{
        for (i = 0; i < DIMENSION; i++){
            axis_ptr->direction[i] = middle[i] / r;
        }
    }

    if (!this->establish_pairs(axis_ptr)){
        return 0;
    }
    if (!this->check_transform_order(axis_ptr)){
        return 0;
    }
    if (!this->optimize_transform_params(axis_ptr)){
        return 0;
    }
    
    return 1;
}



bool Symmetry::detect_rot_sym(std::shared_ptr<Structure> mol, std::vector<int> torsion_atoms, int order) // ANGLE STATT ORDER
{   
    this->support_atom = torsion_atoms[0];
    this->support = mol->atoms[this->support_atom]->coords;
    this->AtomsCount = torsion_atoms.size();

    for (int atom: torsion_atoms){
        this->atoms.push_back(std::make_shared<Atom>(*(mol->atoms[atom])));
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
