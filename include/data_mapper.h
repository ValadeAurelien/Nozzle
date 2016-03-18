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

// ========= obj dependencies ========

class Usr_Interface;
struct arglist_struct;
struct mesh_struct;
typedef vector<vector<mesh_struct>> mesh_grid_t;

// ======== interface ========

class Data_Mapper
{
        public:
                Data_Mapper(Usr_Interface *UI, arglist_struct *arglist_pt);
                void create_datafile_from_mesh_grid(mesh_grid_t *mesh_grid_pt);
                void show_datafile();
        private:
                const string PATH = "./datafiles/";
                Usr_Interface *UI;
                arglist_struct *arglist_pt;
};


#endif


