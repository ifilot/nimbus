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

#include "hf.h"

HF::HF(const std::shared_ptr<Molecule>& _mol) :
mol(_mol) {
    this->cgfs = this->mol->get_cgfs();
}

void HF::scf() {
    // Construct matrices to hold values for the overlap,
    // kinetic and two nuclear integral values, respectively.
    Eigen::MatrixXd S = Eigen::MatrixXd::Zero(this->cgfs->size(), this->cgfs->size());
    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(this->cgfs->size(), this->cgfs->size());
    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(this->cgfs->size(), this->cgfs->size());

    // calculate the integral values using the integrator class
    for(unsigned int i=0; i<this->cgfs->size(); i++) {
        for(unsigned int j=i; j<this->cgfs->size(); j++) {
            S(i,j) = S(j,i) = integrator.overlap(this->cgfs->at(i), this->cgfs->at(j));
            T(i,j) = T(j,i) = integrator.kinetic(this->cgfs->at(i), this->cgfs->at(j));

            for(unsigned int k=0; k<this->mol->get_nr_atoms(); k++) {
                double v = integrator.nuclear(this->cgfs->at(i),
                                              this->cgfs->at(j),
                                              this->mol->get_atomic_position(k),
                                              this->mol->get_atomic_charge(k));
                V(i,j) += v;
            }

            V(j,i) = V(i,j);
        }
    }

    // Calculate 1-electron Hamiltonian matrix
    Eigen::MatrixXd H = T + V;

    // calculate all two-electron integrals
    unsigned int size = integrator.teindex(this->cgfs->size(),this->cgfs->size(),this->cgfs->size(),this->cgfs->size());
    std::vector<double> tedouble(size, -1.0);
    for(unsigned int i=0; i<this->cgfs->size(); i++) {
        for(unsigned int j=0; j<this->cgfs->size(); j++) {
            unsigned int ij = i*(i+1)/2 + j;
            for(unsigned int k=0; k<this->cgfs->size(); k++) {
                for(unsigned int l=0; l<this->cgfs->size(); l++) {
                    unsigned int kl = k * (k+1)/2 + l;
                    if(ij <= kl) {
                        unsigned int idx = integrator.teindex(i,j,k,l);
                        tedouble[idx] = integrator.repulsion(this->cgfs->at(i), this->cgfs->at(j), this->cgfs->at(k), this->cgfs->at(l));
                    }
                }
            }
        }
    }

    // perform a canonical diagonalization on the overlap matrix to
    // obtain orthonormal spinorbitals (required for the Slater
    // Determinant)
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S);
    Eigen::MatrixXd D = es.eigenvalues().real().asDiagonal();
    Eigen::MatrixXd U = es.eigenvectors().real();

    for(unsigned int i=0; i<this->cgfs->size(); i++) {
        D(i,i) = 1.0 / sqrt(D(i,i));
    }


    // Calculate the transformation matrix
    Eigen::MatrixXd X = U * D;
    Eigen::MatrixXd Xp = X.transpose();

    /*********************************
     *
     * STEP 4: Obtain initial guess for density matrix
     *
     *********************************/

    // Create Density matrices
    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(this->cgfs->size(), this->cgfs->size());
    Eigen::MatrixXd Pnew = Eigen::MatrixXd::Zero(this->cgfs->size(), this->cgfs->size());

    // Create two-electron Hamiltonian matrix
    Eigen::MatrixXd G = Eigen::MatrixXd::Zero(this->cgfs->size(), this->cgfs->size());

    // define mixing parameter in the SCF iterations
    static const double alpha = 0.5;

    // hold energy value of previous iteration
    double energy_old = 0.0;

    // difference between previous and current SCF iteration
    // (initialize with some large number)
    double energy_difference = 1.0;

    // keep track of number of iterations
    unsigned int loop_counter = 0;

    // construct object to store eigenvectors and eigenvalues of the Fock-matrix in
    Eigen::VectorXd orbital_energies;

    // calculate nuclear repulsion
    double e_nuc = 0.0;
    for(unsigned int i=0; i<this->mol->get_nr_atoms(); i++) {
        for(unsigned int j=i+1; j<this->mol->get_nr_atoms(); j++) {
            e_nuc += this->mol->get_atomic_charge(i) * this->mol->get_atomic_charge(j) /
                     (this->mol->get_atomic_position(i) - this->mol->get_atomic_position(j)).norm();
        }
    }

    /*
     * START ITERATIVE PROCEDURE
     */
    std::cout << "-------------------------------------------------" << std::endl;
    while(energy_difference > 1e-5 && loop_counter < 100) {
        loop_counter++; // increment loop counter

        /*********************************
         *
         * STEP 5: Calculate G, H, F and F' from P
         *
         *********************************/

        // Populate two-electron hamiltonian matrix
        for(unsigned int i=0; i<this->cgfs->size(); i++) {
            for(unsigned int j=0; j<this->cgfs->size(); j++) {
                G(i,j) = 0.; /* reset G matrix */
                for(unsigned int k=0; k<this->cgfs->size(); k++) {
                    for(unsigned int l=0; l<this->cgfs->size(); l++) {
                        unsigned int index1 = integrator.teindex(i,j,l,k);
                        unsigned int index2 = integrator.teindex(i,k,l,j);
                        G(i,j) += P(k,l) * (tedouble[index1] - 0.5 * tedouble[index2]);
                    }
                }
            }
        }

        // Calculate Fock Matrix
        Eigen::MatrixXd F = H + G;

        // Transform Fock Matrix using our basis transformation matrix
        Eigen::MatrixXd Fp = Xp * F * X;


        /*********************************
         *
         * STEP 6: Diagonalize F' to obtain C' and e
         *
         *********************************/

        // Calculate eigenvalues and vectors
        es.compute(Fp);
        Eigen::MatrixXd Cc = es.eigenvectors().real();
        Eigen::MatrixXd en = es.eigenvalues().real().asDiagonal();

        // Calculate energy
        double energy = 0.0;
        Eigen::MatrixXd M = H + F;
        for(unsigned int i=0; i<this->cgfs->size(); i++) {
            for(unsigned int j=0; j<this->cgfs->size(); j++) {
                energy += P(j,i) * M(i,j);
            }
        }

        // Add nuclear repulsion to the orbital energies
        energy = energy * 0.5 + e_nuc;

        /*********************************
         *
         * STEP 7: Diagonalize C from C'
         *
         *********************************/

        // Obtain true coefficient matrix using the transformation matrix
        this->C = X * Cc;
        orbital_energies = en.diagonal();

        /*********************************
         *
         * STEP 8: Calculate new P from C
         *
         *********************************/

        // obtain new P matrix from the old matrix
        Pnew = Eigen::MatrixXd::Zero(this->cgfs->size(), this->cgfs->size());
        for(unsigned int i=0; i<this->cgfs->size(); i++) {
            for(unsigned int j=0; j<this->cgfs->size(); j++) {
                for(unsigned int k=0; k<this->mol->get_nr_elec() / 2; k++) {
                    Pnew(i,j) += 2.0 * C(i,k) * C(j,k);
                }
            }
        }

        // construct new P-matrix by mixing old and new matrix. This is not always necessary,
        // but sometimes the SCF calculation does not converge without it.
        for(unsigned int i=0; i<this->cgfs->size(); i++) {
            for(unsigned int j=0; j<this->cgfs->size(); j++) {
                P(i,j) = (1.0-alpha) * Pnew(i,j) + alpha * P(i,j);
            }
        }

        // report energy and calculate difference with previous energy value
        energy_difference = std::abs(energy - energy_old);
        energy_old = energy;

        std::cout << (boost::format("%4i | %12.8f") % (loop_counter) % energy).str() << std::endl;
    }

    std::cout << "-------------------------------------------------" << std::endl;
    std::cout << "Stopping because energy convergence was achieved." << std::endl;
    std::cout << std::endl;

    // outputting orbital energies
    std::cout << "-----+--------------" << std::endl;
    std::cout << "  #  | Energy" << std::endl;
    std::cout << "-----+--------------" << std::endl;
    for(unsigned int i=0; i<10; i++) {
        std::cout << (boost::format("%4i | %12.8f") % (i+1) % orbital_energies[i]).str() << std::endl;
    }
    std::cout << std::endl;
}

void HF::write_charge_files() {
    // construct density plotter object and export the electron density of
    // each molecular orbital in our system to a seperate file
    DensityPlotter dp;

    // output the molecular orbital density to a set of files
    dp.plot_densities_chgcar(*this->cgfs, this->C, this->mol->get_nr_elec() / 2 + 1);
}
