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
#include "diff_eq_solver.h"

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
        system( ("mkdir 2>/dev/null " + this->data_dir).c_str() );      
}


void Data_Mapper::create_datafile_from_mesh_grid(Diff_Eq_Solver *DES)
{
        mesh_grid_t &mesh_grid = (*(DES->mesh_grid_pt2));

        this->create_directory();
        this->UI->cout_str("DM-> saving mesh grid datas in " + this->data_dir);

        ofstream file;
        file.open(this->data_dir + "/meshgrid.data", ios::out);

        file << "# Result of the experimentation: " << this->argfile_name << endl << endl;
        file << "# Maximums (asbolute values): " << endl
            << "# pressure = " << setprecision(PRECISION) << DES->variables_max.pressure << endl
            << "# temperature = " << setprecision(PRECISION) << DES->variables_max.temperature << endl
            << "# vol_mass = " << setprecision(PRECISION) << DES->variables_max.vol_mass << endl
            << "# speed0 = " << setprecision(PRECISION) << DES->variables_max.speed0 << endl
            << "# speed1 = " << setprecision(PRECISION) << DES->variables_max.speed1 << endl
            << "# turb_en = " << setprecision(PRECISION) << DES->variables_max.turb_en << endl
            << "# turb_dis = " << setprecision(PRECISION) << DES->variables_max.turb_dis << endl
            << endl;

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

void Data_Mapper::thrust_plotter(Diff_Eq_Solver *DES) {
        
        this->create_directory();
        this->UI->cout_str("DM-> saving time step and thrust datas in " + this->data_dir);

        ofstream file;
        file.open(this->data_dir + "/time_step_n_thrust.data", ios::out);
        file << "# Result of the experimentation: " << this->argfile_name << endl << endl;
        file << "#" << setw(WIDTH) << "time" <<  setw(WIDTH) << "time_step" << setw(WIDTH) << "thrust" << endl << endl;

        double time = 0;
        for (int k=0; k<this->arglist_pt->iter_number_solver;k++) {
                file << setw(WIDTH) << setprecision(PRECISION) << time
                    << setw(WIDTH) << setprecision(PRECISION) << DES->time_steps[k]
                    << setw(WIDTH) << setprecision(PRECISION) << DES->thrust[k] 
                    << endl;
                time += DES->time_steps[k];
        }
        file.close();
}

void show()
{

}
