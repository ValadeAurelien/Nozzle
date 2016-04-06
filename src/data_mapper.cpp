/* This is database mapper implementation linked to data_mapper.h
 *
 *
 *
 *
*/


// ======== includes ======== 
#include "../include/data_mapper.h"
#include "../include/arg_interpreter.h"
#include "../include/usr_interface.h"
#include "../include/nozzle_profiler.h"

#include <iostream>
#include <fstream>
#include <iomanip>

#define WIDTH 18
#define PRECISION 7

// ======== implementation ========

Data_Mapper::Data_Mapper(Usr_Interface *UI, arglist_struct *arglist_pt)
{
        this->UI = UI;
        this->arglist_pt = arglist_pt;
}

void Data_Mapper::create_datafile_from_mesh_grid(mesh_grid_t *mesh_grid_pt)
{
        string datafile_path = this->PATH + this->arglist_pt->datafile_name;
        this->UI->cout_str("DM-> saving mesh grid datas in " + datafile_path);
        
        mesh_grid_t &mesh_grid = *mesh_grid_pt;

        ofstream file;
        file.open(datafile_path, ios::out);
        bool b = true;
        file << "# Result of the experimentation with args: " << b<< endl;
        file << "#" << setw(WIDTH) << "x" << setw(WIDTH) << "y"
                << setw(WIDTH) << "is_wall" << setw(WIDTH) << "pressure"
                << setw(WIDTH) << "temper" << setw(WIDTH) << "vol_mass"
                << setw(WIDTH) << "speed_1" << setw(WIDTH) <<  "speed_2"
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
                                << endl;
                }
                file << endl;
        }
        file.close();
}

void Data_Mapper::thrust_plotter(vector<data_t> *thrust_pt) {
        ofstream file;
        file.open("./thrust.data", ios::out);
        double current_time;
        for (int k=0; k<thrust_pt->size();k++) {
                current_time = k*this->arglist_pt->time_step;
                file << setw(WIDTH) << k*this->arglist_pt->time_step << setw(WIDTH) << setprecision(PRECISION) << (*(thrust_pt))[k] << endl;
        }
        file.close();
}

void show()
{

}
