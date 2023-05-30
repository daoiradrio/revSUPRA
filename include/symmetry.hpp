/*
 *  Brute force symmetry analyzer.
 *  This is actually C++ program, masquerading as a C one!
 *
 *  (C) 1996, 2003 S. Patchkovskii, Serguei.Patchkovskii@sympatico.ca
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 * $Log: symmetry.c,v $
 * Revision 1.16  2003/04/04  13:05:03  patchkov
 * Revision 1.15  2000/01/25  16:47:17  patchkov
 * Revision 1.14  2000/01/25  16:39:08  patchkov
 * Revision 1.13  1996/05/24  12:32:08  ps
 * Revision 1.12  1996/05/23  16:10:47  ps
 * First reasonably stable version.
 *
 * Cloned from GitHub Repository https://github.com/nquesada/symmetry
 */

/*
 *  All specific structures should have corresponding elements in the
 *  same position generic structure does.
 *
 *  Planes are characterized by the surface normal direction
 *  (taken in the direction *from* the coordinate origin)
 *  and distance from the coordinate origin to the plane
 *  in the direction of the surface normal.
 *
 *  Inversion is characterized by location of the inversion center.
 *
 *  Rotation is characterized by a vector (distance+direction) from the origin
 *  to the rotation axis, axis direction and rotation order. Rotations
 *  are in the clockwise direction looking opposite to the direction
 *  of the axis. Note that this definition of the rotation axis
 *  is *not* unique, since an arbitrary multiple of the axis direction
 *  can be added to the position vector without changing actual operation.
 *
 *  Mirror rotation is defined by the same parameters as normal rotation,
 *  but the origin is now unambiguous since it defines the position of the
 *  plane associated with the axis.
 *
 */

#ifndef SYMMETRY_HPP
#define SYMMETRY_HPP

#include <utils.hpp>
#include <structure.hpp>
#include <symmetryelement.hpp>
#include <rotationaxis.hpp>

#include <vector>
#include <memory>
#include <cmath>
#include <math.h>
#include <iostream>
#include <algorithm>

//#ifndef M_PI
//#define M_PI 3.1415926535897932384626433832795028841971694
//#endif

#define DIMENSION 3
#define MAXPARAM  7



class Symmetry{
    public:
        Symmetry();
        ~Symmetry();

        int                                 support_atom;
        int                                 AtomsCount;
        std::vector<double>                 support;
        std::vector<double>                 geom_center;
        std::vector<double>                 dist_geom_center;
        std::vector<std::shared_ptr<Atom>>  atoms;

        bool    detect_rot_sym(std::shared_ptr<Structure> mol, std::vector<int> torsion_atoms, int order);
        void    find_geometric_center();
        void    check_C2_axis();
        int     init_C2(int i, int j);
        int     establish_pairs(std::shared_ptr<SymmetryElement> elem);
        int     check_transform_order(std::shared_ptr<SymmetryElement> elem);
        int     optimize_transform_params(std::shared_ptr<SymmetryElement> elem);
        double  eval_opt_target_func(std::shared_ptr<SymmetryElement> elem, std::shared_ptr<int> finish);
        void    get_params(std::shared_ptr<SymmetryElement> elem, std::vector<double> values);
        void    set_params(std::shared_ptr<SymmetryElement> elem, std::vector<double> values);

    private:

};



#endif
