/* This is datamapper header
 *
 * This module has to help to browse and build the database
 *
 * The object defined below is constructed in main.cpp and called through its adress elsewhere
 *
*/

// ======== include guards ========

#ifndef __DATA_MAPPER_IS_INCLUDED__
#define __DATA_MAPPER_IS_INCLUDED__

using namespace std;

// ======= includes ========

#include <string>
#include <vector>
#include <fstream>

// ========= obj dependencies ========


class Usr_Interface;
class Diff_Eq_Solver;
struct arglist_struct;
struct mesh_struct;
typedef vector<vector<mesh_struct>> mesh_grid_t;
typedef long double data_t;

// ======== interface ========

class Data_Mapper
{
        public:
                Data_Mapper(Usr_Interface *UI, arglist_struct *arglist_pt, string argfile_name);
                void create_directory();
                void write_args();
                void create_datafile_from_mesh_grid(Diff_Eq_Solver * DES);
                void show_datafile();
                void temporal_stuff_plotter(Diff_Eq_Solver * DES);
        private:
                const string PATH = "./datafiles/";
                Usr_Interface *UI;
                arglist_struct *arglist_pt;
                string argfile_name;
                string data_dir;
                ofstream file;
};


#endif


