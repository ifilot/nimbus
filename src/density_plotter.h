/**************************************************************************
 *                                                                        *
 *   Author: Ivo Filot <i.a.w.filot@tue.nl>                               *
 *                                                                        *
 *   NIMBUS is free software:                                             *
 *   you can redistribute it and/or modify it under the terms of the      *
 *   GNU General Public License as published by the Free Software         *
 *   Foundation, either version 3 of the License, or (at your option)     *
 *   any later version.                                                   *
 *                                                                        *
 *   NIMBUS is distributed in the hope that it will be useful,            *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty          *
 *   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.              *
 *   See the GNU General Public License for more details.                 *
 *                                                                        *
 *   You should have received a copy of the GNU General Public License    *
 *   along with this program.  If not, see http://www.gnu.org/licenses/.  *
 *                                                                        *
 **************************************************************************/

#ifndef _DENSITY_PLOTTER_H
#define _DENSITY_PLOTTER_H

#include <fstream>
#include <string>
#include <Eigen/Eigenvalues>
#include <boost/format.hpp>
#include <mutex>

#include "gamma.h"
#include "cgf.h"
#include "molecule.h"

/**
 * @brief      Class that is able to export electron densities of the molecular
 *             orbitals to a series of files for visualization
 */
class DensityPlotter {
private:
    std::shared_ptr<Molecule> mol;

public:
    /**
     * @brief      default constructor
     */
    DensityPlotter(const std::shared_ptr<Molecule>& _mol);

    /**
     * @brief      export electron densities of each MO to a CHGCAR file (VASP
     *             density format)
     *
     * @param[in]  cgfs        reference to vector of contracted gaussian
     *                         functionals
     * @param[in]  C           reference to coefficient matrix C
     * @param[in]  nr          number of orbitals
     * @param[in]  boxsize     size of the unitcell
     * @param[in]  resolution  density object resolution
     */
    void plot_densities_chgcar(const std::vector<CGF>& cgfs,
                               const Eigen::MatrixXd& C,
                               unsigned int nr,
                               double boxsize,
                               double resolution);

    /**
     * @brief      export electron densities of each MO to a cube file (Gaussian
     *             format)
     *
     * @param[in]  atoms  vector of atom number of atoms
     * @param[in]  pos    position of the atoms
     * @param[in]  cgfs   set of contracted gaussian functions (basis set)
     * @param[in]  C      coefficient matrix
     */
    void plot_densities_cube(const std::vector<unsigned int>& atoms,
                             const std::vector<vec3>& pos,
                             const std::vector<CGF>& cgfs,
                             const Eigen::MatrixXd& C);

private:
    /**
     * @brief      generate list of atoms for charge file
     *
     * @return     string to print
     */
    std::string generate_atom_list();

    /**
     * @brief      generate list of atom coordinates for charge file
     *
     * @param[in]  boxsize  size of the unit cell
     *
     * @return     string to print
     */
    std::string generate_atom_coordinates(double boxsize);

};

#endif // _DENSITY_PLOTTER_H
