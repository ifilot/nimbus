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

#include "density_plotter.h"

/**
 * @brief      default constructor
 */
DensityPlotter::DensityPlotter(const std::shared_ptr<Molecule>& _mol) :
mol(_mol) {}

/**
 * @brief      export electron densities of each MO to a CHGCAR file (VASP density format)
 *
 * @param[in]  cgfs  reference to vector of contracted gaussian functionals
 * @param[in]  C     reference to coefficient matrix C
 */
void DensityPlotter::plot_densities_chgcar(const std::vector<CGF>& cgfs,
                                           const Eigen::MatrixXd& C,
                                           unsigned int nr,
                                           double boxsize,
                                           double resolution) {
    const double range = boxsize / 2.0;
    const double inc = resolution;
    const double volume = boxsize * boxsize * boxsize;

    std::cout << "Outputting density files:" << std::endl;
    std::cout << "-------------------------" << std::endl;

    const unsigned int sz = (int)(range/inc*2+1);

    std::vector<std::vector<double> > cache;
    for(unsigned int j=0; j<cgfs.size(); j++) {
        cache.push_back(std::vector<double>(sz * sz *sz));
    }

    std::cout << "Pre-caching wavefunction values...";
    auto start = std::chrono::system_clock::now();

    #pragma omp parallel for
    for(unsigned int zz = 0; zz < sz; zz++) {
        double z = -range + zz * inc;

        for(unsigned int yy = 0; yy < sz; yy++) {
            double y = -range + yy * inc;

            for(unsigned int xx = 0; xx < sz; xx++) {
                double x = -range + xx * inc;

                unsigned int idx = zz * sz * sz + yy * sz + xx;

                for(unsigned int j=0; j<cgfs.size(); j++) {
                    cache[j][idx] = cgfs[j].get_value(x,y,z);
                }
            }
        }
    }

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "[" << elapsed_seconds.count() << " seconds]" << std::endl;

    std::cout << "Using boxsize of " << boxsize << " A and resolution of " << resolution << " A." << std::endl;

    // calculate some listings
    const std::string atomlist = this->generate_atom_list();
    const std::string coordinatelist = this->generate_atom_coordinates(boxsize);

    for(unsigned int i=0; i<nr; i++) {
        Eigen::VectorXd coefs = C.col(i);

        std::string filename = (boost::format("%04i.dat") % (i+1)).str();
        std::cout << "Writing " << filename << std::endl;
        std::ofstream outfile(filename);

        outfile << "Amplitude of Molecular Orbital #" << (i+1) << std::endl;
        outfile << "1.0" << std::endl;
        outfile << boost::format("%12.6f  %12.6f  %12.6f\n") % boxsize % 0.0 % 0.0;
        outfile << boost::format("%12.6f  %12.6f  %12.6f\n") % 0.0 % boxsize % 0.0;
        outfile << boost::format("%12.6f  %12.6f  %12.6f\n") % 0.0 % 0.0 % boxsize;
        outfile << atomlist;
        outfile << coordinatelist;
        outfile << "" << std::endl;
        outfile << sz << " ";
        outfile << sz << " ";
        outfile << sz << std::endl;

        std::vector<double> values(sz * sz * sz, 0.0);
        unsigned int cnt = 0;

        #pragma omp parallel for
        for(unsigned int zz = 0; zz < sz; zz++) {
            double z = -range + zz * inc;

            for(unsigned int yy = 0; yy < sz; yy++) {
                double y = -range + yy * inc;

                for(unsigned int xx = 0; xx < sz; xx++) {
                    double x = -range + xx * inc;

                    unsigned int idx = zz * sz * sz + yy * sz + xx;

                    for(unsigned int j=0; j<cgfs.size(); j++) {
                        values[idx] += coefs[j] * cache[j][idx];
                    }
                }
            }
        }

        size_t nrthreads = omp_get_max_threads();
        omp_set_num_threads(nrthreads); // always allocate max threads
        std::stringstream local[nrthreads];

        #pragma omp parallel
        {
            size_t threadnum = omp_get_thread_num();

            // calculate size
            size_t data = values.size() / nrthreads;
            size_t rem = values.size() % nrthreads;

            // divide task
            size_t start = values.size() / nrthreads * threadnum;
            size_t stop = values.size() / nrthreads * (threadnum + 1);
            if(threadnum == nrthreads - 1) {
                stop += rem;
            }

            char buffer[50];
            unsigned int cnt = 0;

            for(size_t i=start; i<stop; i++) {

                sprintf(buffer, "% 11.10E", values[i] * volume);
                local[threadnum] << buffer;

                if((i+1) % 5 == 0) {
                    local[threadnum] << "\n";
                } else {
                    local[threadnum] << " ";
                }

                if(threadnum == 0) {
                    if(cnt == 1000) {
                        cnt = 0;
                    }
                    cnt++;
                }
            }
        }

        for(unsigned int i=0; i<nrthreads; i++) {
            outfile << local[i].str();
        }

        outfile.close();
    }
}

/**
 * @brief      export electron densities of each MO to a cube file (Gaussian
 *             format)
 *
 * @param[in]  atoms  vector of atom number of atoms
 * @param[in]  pos    position of the atoms
 * @param[in]  cgfs   set of contracted gaussian functions (basis set)
 * @param[in]  C      coefficient matrix
 */
void DensityPlotter::plot_densities_cube(const std::vector<unsigned int>& atoms,
                                         const std::vector<vec3>& pos,
                                         const std::vector<CGF>& cgfs,
                                         const Eigen::MatrixXd& C) {
    static const double range = 10;
    static const double inc = 0.2;

    for(unsigned int i=0; i<cgfs.size(); i++) {
        Eigen::VectorXd coefs = C.col(i);

        std::string filename = (boost::format("%04i.cube") % (i+1)).str();
        std::cout << "Outputting cube file to: " << filename << std::endl;
        std::ofstream outfile(filename);

        outfile << "HFHSL CUBE FILE." << std::endl;
        outfile << "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z" << std::endl;
        // print number of atoms and center of unit cell
        outfile << (boost::format("%i  %g  %g  %g") % atoms.size() % range % range % range).str() << std::endl;

        // print unit cell information
        outfile << (boost::format("%i  %g  %g  %g") % (int)(range/inc*2+1) % inc % 0.0 % 0.0).str() << std::endl;
        outfile << (boost::format("%i  %g  %g  %g") % (int)(range/inc*2+1) % 0.0 % inc % 0.0).str() << std::endl;
        outfile << (boost::format("%i  %g  %g  %g") % (int)(range/inc*2+1) % 0.0 % 0.0 % inc).str() << std::endl;

        // print atom positions
        for(unsigned int i=0; i<atoms.size(); i++) {
            outfile << (boost::format("%i  %g  %g  %g") % atoms[i] % pos[i][0] % pos[i][1] % pos[i][2]).str() << std::endl;
        }

        unsigned int cnt = 0;
        for(double x=-range; x<=range; x+=inc) {
            for(double y=-range; y<=range; y+=inc) {
                for(double z=-range; z<=range; z+=inc) {

                    double sum = 0.0;
                    for(unsigned int j=0; j<cgfs.size(); j++) {
                        sum += coefs[j] * cgfs[j].get_value(x,y,z);
                    }

                    outfile << (boost::format("  %g") % ((sum * sum) * 2.0)).str();

                    cnt++;
                    if(cnt % 6 == 0) {
                        outfile << std::endl;
                        cnt = 0;
                    }
                }
                outfile << std::endl;
                cnt = 0;
            }
        }
    }
}

/**
 * @brief      generate list of atoms for charge file
 *
 * @return     string to print
 */
std::string DensityPlotter::generate_atom_list() {
    static const std::vector<std::string> names = {"H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na",
                                                   "Mg", "Al", "Si", "P", "S", "Cl", "Ar"};

    std::vector<unsigned int> elnrs;
    std::vector<unsigned int> nrs;
    for(unsigned int i=0; i<this->mol->get_nr_atoms(); i++) {
        unsigned int elnr = this->mol->get_atom(i)->get_charge();
        bool flag = false;

        for(unsigned int j=0; j<elnrs.size(); j++) {
            if(elnrs[j] == elnr) {
                flag = true;
                nrs[j]++;
                break;
            }
        }

        if(!flag) {
            elnrs.push_back(elnr);
            nrs.push_back(1);
        }
    }

    std::stringstream str;
    for(unsigned int i=0; i<elnrs.size(); i++) {
        str << "  " << names[elnrs[i]-1];
    }
    str << std::endl;

    for(unsigned int i=0; i<nrs.size(); i++) {
        str << "  " << nrs[i];
    }
    str << std::endl;

    return str.str();
}

/**
 * @brief      generate list of atom coordinates for charge file
 *
 * @param[in]  boxsize  size of the unit cell
 *
 * @return     string to print
 */
std::string DensityPlotter::generate_atom_coordinates(double boxsize) {
    std::vector<unsigned int> elnrs;
    std::vector<unsigned int> nrs;
    for(unsigned int i=0; i<this->mol->get_nr_atoms(); i++) {
        unsigned int elnr = this->mol->get_atom(i)->get_charge();
        bool flag = false;

        for(unsigned int j=0; j<elnrs.size(); j++) {
            if(elnrs[j] == elnr) {
                flag = true;
                nrs[j]++;
                break;
            }
        }

        if(!flag) {
            elnrs.push_back(elnr);
            nrs.push_back(1);
        }
    }

    std::stringstream str;
    str << "Direct" << std::endl;

    for(unsigned int i=0; i<elnrs.size(); i++) {
        for(unsigned int j=0; j<this->mol->get_nr_atoms(); j++) {
            if(elnrs[i] == this->mol->get_atom(j)->get_charge()) {
                vec3 pos = this->mol->get_atom(j)->get_position() / boxsize;
                str << boost::format("  %12.8f  %12.8f  %12.8f\n") % (pos[0]+0.5) % (pos[1]+0.5) % (pos[2]+0.5);
            }
        }
    }

    return str.str();
}
