/* This is database mapper implementation linked to data_mapper.h
 *
 *
 *
 *
*/


// ======== includes ======== 
#include "data_mapper.h"
#include "arg_interpreter.h"
#include "usr_interface.h"
#include "nozzle_profiler.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>

#define WIDTH 18
#define PRECISION 7

// ======== implementation ========

Data_Mapper::Data_Mapper(Usr_Interface *UI, arglist_struct *arglist_pt, string argfile_name)
{
        this->UI = UI;
        this->arglist_pt = arglist_pt;
        this->argfile_name = argfile_name;
        string::size_type pos = argfile_name.find(".arg");
        string sub_data_dir = argfile_name.substr(0, pos);

        this->data_dir = this->PATH + sub_data_dir;

}

void Data_Mapper::create_directory()
{
        this->UI->cout_str("DM-> saving mesh grid datas in " + this->data_dir);
        system( ("mkdir 2>/dev/null " + this->data_dir).c_str() );      
}


void Data_Mapper::create_datafile_from_mesh_grid(mesh_grid_t *mesh_grid_pt)
{
        mesh_grid_t &mesh_grid = *mesh_grid_pt;

        this->create_directory();

        ofstream file;
        file.open(this->data_dir + "/final_prof_meshgrid.data", ios::out);

        file << "# Result of the experimentation: " << this->argfile_name << endl << endl;
        file << "#" << setw(WIDTH) << "x" << setw(WIDTH) << "y"
                << setw(WIDTH) << "is_wall" << setw(WIDTH) << "pressure"
                << setw(WIDTH) << "temper" << setw(WIDTH) << "vol_mass"
                << setw(WIDTH) << "speed_1" << setw(WIDTH) <<  "speed_2"
                << setw(WIDTH) << "turb_en" << setw(WIDTH) << "turb_dis"
                << endl << endl;

        register int i,j;
        for (i=0; i<this->arglist_pt->x_size; i++){
                for (j=0; j<this->arglist_pt->y_size; j++){
                        file << setw(WIDTH) << setprecision(PRECISION) << i*this->arglist_pt->space_step
                                << setw(WIDTH) << setprecision(PRECISION) << j*this->arglist_pt->space_step
                                << setw(WIDTH) << setprecision(PRECISION) << mesh_grid[i][j].is_wall
                                << setw(WIDTH) << setprecision(PRECISION) << mesh_grid[i][j].pressure
                                << setw(WIDTH) << setprecision(PRECISION) << mesh_grid[i][j].temperature
                                << setw(WIDTH) << setprecision(PRECISION) << mesh_grid[i][j].vol_mass
                                << setw(WIDTH) << setprecision(PRECISION) << mesh_grid[i][j].speed[0]
                                << setw(WIDTH) << setprecision(PRECISION) << mesh_grid[i][j].speed[1]
                                << setw(WIDTH) << setprecision(PRECISION) << mesh_grid[i][j].turb_en
                                << setw(WIDTH) << setprecision(PRECISION) << mesh_grid[i][j].turb_dis
                                << endl;
                }
                file << endl;
        }
        file.close();
}

void Data_Mapper::thrust_plotter(vector<data_t> *thrust_pt) {
        
        this->create_directory();
        
        ofstream file;
        file.open(this->data_dir + "/final_prof_thrust.data", ios::out);
        for (int k=0; k<thrust_pt->size();k++) {
                file << setw(WIDTH) << k*this->arglist_pt->time_step 
                    << setw(WIDTH) << setprecision(PRECISION) << (*(thrust_pt))[k] 
                    << endl;
        }
        file.close();
}

void show()
{

}
