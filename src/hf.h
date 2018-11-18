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

#ifndef _HF_H
#define _HF_H

#include <Eigen/Eigenvalues>
#include <memory>
#include <boost/format.hpp>

#include "integrals.h"
#include "density_plotter.h"
#include "cgf.h"
#include "molecule.h"

class HF {
private:
    const std::vector<CGF>* cgfs;
    Integrator integrator;
    std::shared_ptr<Molecule> mol;
    Eigen::MatrixXd S;
    Eigen::MatrixXd C;
    Eigen::VectorXd orbital_energies;

public:
    HF(const std::shared_ptr<Molecule>& _mol);

    /**
     * @brief      perform self-consistent field calculation
     */
    void scf();

    /**
     * @brief      Writes wavefunction amplitude files
     *
     * @param[in]  boxsize     boxsize
     * @param[in]  resolution  resolution
     */
    void write_charge_files(double boxsize, double resolution);

    /**
     * @brief      Write coefficient matrix and orbital energies to file
     *
     * @param[in]  filename  The filename
     */
    void write_matrices(const std::string& filename);

private:

};

#endif // _HF_H
