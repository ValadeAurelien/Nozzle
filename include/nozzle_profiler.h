/* This is nozzle profiler header
 *
 * This module is the one which choose the nozzle profile and launch the experimentations
 *
 * The object defined below is constructed in main.cpp and is called through its adress elsewhere
 *
*/

// ======== include guards ========

#ifndef __NOZZLE_PROFILER_IS_INCLUDED__
#define __NOZZLE_PROFILER_IS_INCLUDED__

// ========= obj dependencies ========

class Usr_Interface;
class Data_Mapper;
struct arglist_struct;

// ======= includes ========
#include "diff_eq_solver.h"
#include <vector>


using namespace std;
// ======== structs ========


typedef double data_t;

// template <typename data_t>; //creation du template mesh sur la base d'un data_type float, double, long double
struct mesh_struct{
        bool is_wall;
        data_t pressure;
        data_t temperature;
        data_t vol_mass;
        data_t speed[2];
        data_t turb_en;
        data_t turb_dis;
};

// template <typename data_t>;  // pas sur que cette ligne soit utile
//using mesh_grid_t = vector<vector<mesh<data_t>>; // tableau des mailles de type data_t
typedef vector<vector<mesh_struct>> mesh_grid_t;

// ======== interface ========

// template <typename data_t>;
class Nozzle_Profiler
{
        public:
                //fonctions appelées depuis main

                Nozzle_Profiler(Usr_Interface *UI, Data_Mapper *DM, arglist_struct *arglist_pt); //Constructeur

                void profile(); //switch sur l'alfo de profilage a choisir 
                                

                //fonctions internes 

                void create_mesh_grids();
                void set_init_conditions();
                void save_mesh_grid();

                void set_wall(int i,int j);
                bool is_in_x_range(int i);
        
                // fonction segment
                void init_profile_segment();
                void one_iteration_segment();

                // debug

                void constante();
        private:
                Usr_Interface *UI;
                Data_Mapper *DM;
                Diff_Eq_Solver *DES;
                arglist_struct *arglist_pt;
                mesh_grid_t mesh_grid_1;
                mesh_grid_t mesh_grid_2;
                //mesh_grid_t<data_t> *mesh_grid;
};


#endif


