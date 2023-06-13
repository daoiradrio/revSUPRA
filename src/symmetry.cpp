#include <symmetry.hpp>



Symmetry::Symmetry()
{
    this->possible_orders = {
        6, // corresponds to angle increment  60°
        4, // corresponds to angle increment  90°
        3, // corresponds to angle increment 120°
        2, // corresponds to angle increment 180°
    };
}



Symmetry::~Symmetry(){}



/*int Symmetry::old_establish_pairs(std::shared_ptr<SymmetryElement> elem)
{
    int                     i, j, k;
    int                     best_j;
    double                  distance, best_distance;
    std::vector<int>        atom_used(this->AtomsCount, 0);
    std::shared_ptr<Atom>   symmetric;

    for (i = 0; i < this->AtomsCount; i++){
        //std::cout << i << std::endl;
        if (elem->transform[i] >= this->AtomsCount){
            symmetric = elem->transform_atom(this->atoms[i]);
            //std::cout << symmetric->coords[0] << " " << symmetric->coords[1] << " " << symmetric->coords[2] << std::endl;
            best_j = i;
            best_distance = 2 * TolerancePrimary;
            for (j = 0; j < this->AtomsCount; j++){
                if (this->atoms[i]->pse_num != symmetric->pse_num || atom_used[j]){
                    std::cout << "Continue" << std::endl;
                    continue;
                }
                distance = 0.0;
                for (k = 0; k < DIMENSION; k++){
                    distance += pow(symmetric->coords[k] - this->atoms[j]->coords[k], 2);
                }
                distance = sqrt(distance);
                if (distance < best_distance){
                    best_j = j;
                    best_distance = distance;
                }
            }
            if (best_distance > TolerancePrimary){
                std::cout << "Fail" << std::endl;
                return 0;
            }
            elem->transform[i] = j;
            atom_used[best_j] = 1;
            std::cout << "Geschafft" << std::endl;
        }
    }
    
    return 1;
}



int Symmetry::find_pairs(){
    int                     i, j, k;
    int                     best_j;
    double                  distance, best_distance;
    std::vector<int>        atom_used(this->AtomsCount, 0);
    std::shared_ptr<Atom>   symmetric;

    for (i = 0; i < this->AtomsCount; i++){
        if (this->rot_axis->transform[i] >= this->AtomsCount){
            symmetric = this->rot_axis->transform_atom(this->atoms[i]);
            best_j = i;
            best_distance = 2 * TolerancePrimary;
            for (j = 0; j < this->AtomsCount; j++){
                if (this->atoms[j]->pse_num != symmetric->pse_num || atom_used[j]){
                    continue;
                }
                distance = 0.0;
                for (k = 0; k < DIMENSION; k++){
                    distance += pow(symmetric->coords[k] - this->atoms[j]->coords[k], 2);
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
            this->rot_axis->transform[i] = j;
            atom_used[best_j] = 1;
        }
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
    int                     i;
    int                     cycle = 0;
    int                     hits = 0;
    int                     finish = 0;
    int                     vars = elem->nparam;
    double                  f = 0.0;
    double                  fold, fnew, fnew2, fdn, fup;
    double                  a, b, x;
    double                  snorm;
    std::vector<double>     values;
    std::vector<double>     grad;
    std::vector<double>     force;
    std::vector<double>     step;

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
        this->get_params(elem, values);
        for (i = 0; i < vars; i++){
            values[i] -= GradientStep;
            this->set_params(elem, values);
            fdn = this->eval_opt_target_func(elem, nullptr);
            values[i] += 2*GradientStep;
            this->set_params(elem, values);
            fup = this->eval_opt_target_func(elem, nullptr);
            values[i] -= GradientStep;
            grad[i] = (fup - fdn) / (2*GradientStep);
            force[i] = (fup + fdn - 2*f) / (GradientStep*GradientStep);
        }
        snorm = 0.0;
        for (i = 0; i < vars; i++){
            if (force[i] < 0){
                force[i] = -force[i];
            }
            if (force[i] < 1e-3){
                force[i] = 1e-3;
            }
            else if (force[i] > 1e3){
                force[i] = 1e3;
            }
            step[i] = -grad[i] / force[i];
            snorm += step[i]*step[i];
        }
        snorm = sqrt(snorm);
        if (snorm > MaxOptStep){
            for (i = 0; i < vars; i++){
                step[i] *= MaxOptStep / snorm;
            }
            snorm = MaxOptStep;
        }
        do{
            for (i = 0; i < vars; i++){
                values[i] += step[i];
            }
            this->set_params(elem, values);
            fnew = this->eval_opt_target_func(elem, nullptr);
            if (fnew < f){
                break;
            }
            for (i = 0; i < vars; i++){
                values[i] -= step[i];
                step[i] /= 2;
            }
            this->set_params(elem, values);
            snorm /= 2;
        }while (snorm > MinOptStep);
        if ((snorm > MinOptStep) && (snorm < MaxOptStep / 2)){
            for (i = 0; i < vars; i++){
                values[i] += step[i];
            }
            this->set_params(elem, values);
            fnew2 = this->eval_opt_target_func(elem, nullptr);
            for (i = 0; i < vars; i++){
                values[i] -= 2*step[i];
            }
            a = (4*f - fnew2 - 3*fnew) / 2;
            b = (f + fnew2 - 2*fnew) / 2;
            if (b > 0){
                x = -a/(2*b);
                if (x > 0.2 && x < 1.8){
                    for (i = 0; i < vars; i++){
                        values[i] += x*step[i];
                    }
                }
                else{
                    b = 0;
                }
            }
            if (b <= 0){
                if (fnew2 < fnew){
                    for (i = 0; i < vars; i++){
                        values[i] += 2*step[i];
                    }
                }
                else{
                    for (i = 0; i < vars; i++){
                        values[i] += step[i];
                    }
                }
            }
            this->set_params(elem, values);
        } 
    }while(snorm > MinOptStep && cycle++ < MaxOptCycles);

    f = this->eval_opt_target_func(elem, nullptr);
    if (cycle >= MaxOptCycles){
        return 0;
    }
    return 1;
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
        symmetric = elem->transform_atom(this->atoms[i]);
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
    int i;

    if (values.size() != elem->nparam){
        values.resize(elem->nparam);
    }

    for (i = 0; i < elem->nparam; i++){
        values[i] = elem->distance;
    }

    return;
}



void Symmetry::set_params(std::shared_ptr<SymmetryElement> elem, std::vector<double> values){
    if (values.size() != 0){
        elem->distance = values[0];
    }

    return;
}



void Symmetry::check_C2_axis()
{
    int                 i, j, k;
    int                 sum;
    double              distance_i, distance_j;
    std::vector<int>    flags(this->AtomsCount, 0);

    for (i = 1; i < this->AtomsCount; i++){
        for (j = 0; j < i; j++){
            if (this->atoms[i]->pse_num != this->atoms[j]->pse_num){
                continue;
            }
            if (fabs(this->dist_geom_center[i] - this->dist_geom_center[j]) > TolerancePrimary){
                continue;
            }
            if (this->old_init_C2(i, j)){
                flags[i] = 1;
                flags[j] = 1;
            }
        }
    }

    return;
}



// MUSS AXIS TATSÄCHLICH ALS POINTER VERWENDET WERDEN??
int Symmetry::old_init_C2(int atom1, int atom2)
{   
    std::shared_ptr<RotationAxis>       axis_ptr = std::make_shared<RotationAxis>();
    int                                 i;
    double                              r, r1, r2;
    std::vector<double>                 middle(DIMENSION, 0.0);

    axis_ptr->order = 2;
    axis_ptr->transform.resize(this->AtomsCount);
    std::fill(axis_ptr->transform.begin(), axis_ptr->transform.end(), this->AtomsCount+1);

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

    if (!this->old_establish_pairs(axis_ptr)){
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



int Symmetry::init_rot_axis(int from, int to, int order){
    int                             i;
    double                          r, r1, r2;

    this->rot_axis = std::make_shared<RotationAxis>();
    this->rot_axis->order = order;
    this->rot_axis->axis_from = std::make_shared<Atom>(*(this->atoms[from]));
    this->rot_axis->axis_to = std::make_shared<Atom>(*(this->atoms[to]));
    this->rot_axis->transform.resize(this->AtomsCount);
    std::fill(this->rot_axis->transform.begin(), this->rot_axis->transform.end(), this->AtomsCount+1);

    r = 0.0;

    for (i = 0; i < DIMENSION; i++){
        r += pow(this->geom_center[i], 2);
    }
    r = sqrt(r);
    if (r > 0.0){
        for (i = 0; i < DIMENSION; i++){
            this->rot_axis->normal[i] = this->geom_center[i] / r;
        }
    }
    else{
        this->rot_axis->normal[0] = 1.0;
        this->rot_axis->normal[1] = 0.0;
        this->rot_axis->normal[2] = 0.0;
    }

    this->rot_axis->distance = r;
    r = 0.0;

    for (i = 0; i < DIMENSION; i++){
        this->rot_axis->direction[i] = this->atoms[to]->coords[i] - this->atoms[from]->coords[i];
        r += pow(this->rot_axis->direction[i], 2);
    }
    r = sqrt(r);
    for (i = 0; i < DIMENSION; i++){
        this->rot_axis->direction[i] /= r;
    }
    
    return 0;
}*/



void Symmetry::find_geometric_center()
{
    int                 i, j;
    double              r;
    std::vector<double> coord_sum(DIMENSION, 0.0);

    this->geom_center.resize(DIMENSION);
    this->dist_geom_center.resize(this->n_atoms);

    for (i = 0; i < this->n_atoms; i++){
        for (j = 0; j < DIMENSION; j++){
            coord_sum[j] += this->atoms[i]->coords[j];
        }
    }

    for (i = 0; i < DIMENSION; i++){
        this->geom_center[i] = coord_sum[i] / this->n_atoms;
    }

    for (i = 0; i < this->n_atoms; i++){
        r = 0.0;
        for (j = 0; j < DIMENSION; j++){
            r += pow(this->atoms[i]->coords[j] - this->geom_center[j], 2);
        }
        this->dist_geom_center[i] = r;
    }

    return;
}



bool Symmetry::rot_sym_along_bond(
    const std::shared_ptr<Structure>&   mol,
    const std::vector<int>&             rot_atoms,
    const int&                          axis_from,
    const int&                          axis_to,
    const int&                          order
){   
    int                                 i, j, k;
    int                                 best_j;
    int                                 n_atoms = rot_atoms.size();
    double                              distance, best_distance;
    std::shared_ptr<Atom>               symmetric;
    RotationAxis                        rot_axis(mol->atoms[axis_from], mol->atoms[axis_to]);
    std::vector<int>                    atom_used(rot_atoms.size(), 0);
    std::vector<int>                    transform_pairs(rot_atoms.size(), rot_atoms.size()+1);
    std::vector<std::shared_ptr<Atom>>  atoms;

    for (int atom_i: rot_atoms){
        atoms.push_back(std::make_shared<Atom>(*(mol->atoms[atom_i])));
    }

    for (i = 0; i < n_atoms; i++){
        if (transform_pairs[i] >= n_atoms){
            //std::cout << std::endl;
            //std::cout << i << ": " << atoms[i]->coords[0] << " " << atoms[i]->coords[1] << " " << atoms[i]->coords[2] << std::endl;
            best_j = i;
            best_distance = 1.0;
            //best_distance = 2 * TolerancePrimary;
            symmetric = rot_axis.rotate_atom(atoms[i], 360/order);
            //std::cout << "symmetric: " << symmetric->coords[0] << " " << symmetric->coords[1] << " " << symmetric->coords[2] << std::endl;
            for (j = 0; j < n_atoms; j++){
                if (atoms[j]->pse_num != symmetric->pse_num || atom_used[j]){
                    continue;
                }
                //std::cout << j << ": " << atoms[j]->coords[0] << " " << atoms[j]->coords[1] << " " << atoms[j]->coords[2] << std::endl;
                distance = 0.0;
                for (k = 0; k < DIMENSION; k++){
                    distance += pow(atoms[j]->coords[k] - symmetric->coords[k], 2);
                }
                distance = sqrt(distance);
                //std::cout << "distance: " << distance << std::endl;
                if (distance < best_distance){
                    best_j = j;
                    best_distance = distance;
                }
            }
            if (best_distance > 0.5){
            //if (best_distance > TolerancePrimary){
                //std::cout << "Keine Symmetrie" << std::endl;
                return false;
            }
            transform_pairs[i] = j;
            atom_used[best_j] = 1;
        }
    }

    //std::cout << "Symmetrie" << std::endl;
    return true;
}



int Symmetry::rot_sym_along_bond(
    const std::shared_ptr<Structure>&   mol,
    const std::vector<int>&             rot_atoms, 
    const int&                          from,
    const int&                          to
){
    for (int order: this->possible_orders){
        if (this->rot_sym_along_bond(mol, rot_atoms, from, to, order)){
            return order;
        }
    }

    return 1;
}