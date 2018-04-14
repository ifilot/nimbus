/**************************************************************************
 *                                                                        *
 *   Author: Ivo Filot <i.a.w.filot@tue.nl>                               *
 *                                                                        *
 *   DEN2OBJ is free software:                                            *
 *   you can redistribute it and/or modify it under the terms of the      *
 *   GNU General Public License as published by the Free Software         *
 *   Foundation, either version 3 of the License, or (at your option)     *
 *   any later version.                                                   *
 *                                                                        *
 *   DEN2OBJ is distributed in the hope that it will be useful,           *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty          *
 *   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.              *
 *   See the GNU General Public License for more details.                 *
 *                                                                        *
 *   You should have received a copy of the GNU General Public License    *
 *   along with this program.  If not, see http://www.gnu.org/licenses/.  *
 *                                                                        *
 **************************************************************************/

#include <chrono>
#include <tclap/CmdLine.h>

#include "config.h"
#include "hf.h"

int main(int argc, char* argv[]) {
    try {
        TCLAP::CmdLine cmd("Performs Hartree-Fock simulation on molecule and store all orbital information.", ' ', PROGRAM_VERSION);

        //**************************************
        // declare values to be parsed
        //**************************************

        // input filename
        TCLAP::ValueArg<std::string> arg_input_filename("i","input","Input file (i.e. CHGCAR)",true,"CHGCAR","filename");
        cmd.add(arg_input_filename);

        cmd.parse(argc, argv);

        //**************************************
        // Inform user about execution
        //**************************************
        std::cout << "--------------------------------------------------------------" << std::endl;
        std::cout << "Executing "<< PROGRAM_NAME << " v." << PROGRAM_VERSION << std::endl;
        std::cout << "Author: Ivo Filot <i.a.w.filot@tue.nl>" << std::endl;
        std::cout << "--------------------------------------------------------------" << std::endl;
        std::cout << std::endl;

        //**************************************
        // parsing values
        //**************************************
        std::string input_filename = arg_input_filename.getValue();

        auto start = std::chrono::system_clock::now();

        auto mol = std::make_shared<Molecule>(input_filename);
        HF hf(mol);
        hf.scf();
        hf.write_charge_files();

        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "--------------------------------------------------------------" << std::endl;
        std::cout << "Done in " << elapsed_seconds.count() << " seconds." << std::endl;

        return 0;

    }  catch (TCLAP::ArgException &e) {
        std::cerr << "error: " << e.error() <<
                     " for arg " << e.argId() << std::endl;
        return -1;
    }
}
