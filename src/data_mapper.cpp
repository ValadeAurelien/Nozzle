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

void Data_Mapper::write_args()
{
        this->file << "# Result of the experimentation: " << this->argfile_name << endl << endl;
        this->file << "# With args:" << endl
            << "# --> x_size = " << this->arglist_pt->x_size << endl
            << "# --> y_size = " << this->arglist_pt->y_size << endl
            << "# --> init_time_step = " << this->arglist_pt->init_time_step << endl
            << "# --> space_step = " << this->arglist_pt->space_step << endl
            << "# --> CBL_cond = " <<  this->arglist_pt->CFL_cond << endl
            << "# --> iter_number_solver = " <<  this->arglist_pt->iter_number_solver << endl
            << "# --> iter_number_profiler = " <<  this->arglist_pt->iter_number_profiler << endl
            << "# --> nb_of_threads = " <<  this->arglist_pt->nb_of_threads << endl
            << "# --> ->nb_pts = " <<  this->arglist_pt->nozzle_fitting_init_arg.nb_pts << endl;
        
        register int i;
        
        for (i=0; i<this->arglist_pt->nozzle_fitting_init_arg.nb_pts; i++){
            this->file << "# --> ->abs = " << this->arglist_pt->nozzle_fitting_init_arg.abscisses[i] << endl;
        }

        for (i=0; i<this->arglist_pt->nozzle_fitting_init_arg.nb_pts; i++){
            this->file << "# --> ->ord = " << this->arglist_pt->nozzle_fitting_init_arg.ordinates[i] << endl;
        }
        
        this->file << "# --> init_cond_type = " << this->arglist_pt->init_cond_type << endl
             << "# --> iter_number_evol_chamber = " << this->arglist_pt->iter_number_evol_chamber << endl
             << "# --> differential_equation_solver_algo = " <<  this->arglist_pt->diff_eq_solver_algo << endl
             << "# --> VDW_a_coef = " <<  this->arglist_pt->VDW_a_coef << endl
             << "# --> VDW_b_coef = " <<  this->arglist_pt->VDW_b_coef << endl
             << "# --> thermal_conduction = " <<  this->arglist_pt->thermal_conduction << endl
             << "# --> lambda = " <<  this->arglist_pt->lambda << endl
             << "# --> mol_mass = " <<  this->arglist_pt->mol_mass << endl
             << "# --> dyn_visc = " <<  this->arglist_pt->dyn_visc << endl
             << "# --> init_chamber_pressure = " <<  this->arglist_pt->init_cond.chamber_pressure << endl
             << "# --> init_chamber_temperature = " << this->arglist_pt->init_cond.chamber_temp << endl
             << "# --> init_chamber_speed = " <<  this->arglist_pt->init_cond.chamber_speed << endl
             << "# --> init_chamber_turb_en = " <<  this->arglist_pt->init_cond.chamber_turb_en << endl
             << "# --> init_chamber_turb_dis = " <<  this->arglist_pt->init_cond.chamber_turb_dis << endl
             << "# --> init_atmosphere_pressure = " <<  this->arglist_pt->init_cond.atmosphere_pressure << endl
             << "# --> init_atmosphere_temperature = " <<  this->arglist_pt->init_cond.atmosphere_temp << endl
             << "# --> init_atmosphere_speed = " <<  this->arglist_pt->init_cond.atmosphere_speed << endl
             << "# --> init_atmosphere_turb_en = " <<  this->arglist_pt->init_cond.atmosphere_turb_en << endl
             << "# --> init_atmosphere_turb_dis = " <<  this->arglist_pt->init_cond.atmosphere_turb_dis << endl 
             << endl;


}

void Data_Mapper::create_datafile_from_mesh_grid(Diff_Eq_Solver *DES)
{

        mesh_grid_t &mesh_grid = (*(DES->mesh_grid_pt2));

        this->create_directory();
        this->UI->cout_str("DM-> saving mesh grid datas in " + this->data_dir);

        this->file.open(this->data_dir + "/meshgrid_" + to_string(DES->ite_count + 1) +"_ite.data", ios::out);

        this->write_args();
  
        this->file << "# Maximums (asbolute values) after " << DES->ite_count + 1 << " iterations: " << endl
            << "# pressure = " << setprecision(PRECISION) << DES->variables_max.pressure[DES->ite_count] << endl
            << "# temperature = " << setprecision(PRECISION) << DES->variables_max.temperature[DES->ite_count] << endl
            << "# vol_mass = " << setprecision(PRECISION) << DES->variables_max.vol_mass[DES->ite_count] << endl
            << "# speed0 = " << setprecision(PRECISION) << DES->variables_max.speed0[DES->ite_count] << endl
            << "# speed1 = " << setprecision(PRECISION) << DES->variables_max.speed1[DES->ite_count] << endl
            << "# turb_en = " << setprecision(PRECISION) << DES->variables_max.turb_en[DES->ite_count] << endl
            << "# turb_dis = " << setprecision(PRECISION) << DES->variables_max.turb_dis[DES->ite_count] << endl
            << endl;

        this->file << "#" << setw(WIDTH) << "x" << setw(WIDTH) << "y"
                << setw(WIDTH) << "is_wall" << setw(WIDTH) << "pressure"
                << setw(WIDTH) << "temper" << setw(WIDTH) << "vol_mass"
                << setw(WIDTH) << "speed_1" << setw(WIDTH) <<  "speed_2"
                << setw(WIDTH) << "turb_en" << setw(WIDTH) << "turb_dis"
                << endl << endl;

        register int i,j;
        for (i=0; i<this->arglist_pt->x_size; i++){
                for (j=0; j<this->arglist_pt->y_size; j++){
                        this->file << setw(WIDTH) << setprecision(PRECISION) << i*this->arglist_pt->space_step
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
                this->file << endl;
        }
        this->file.close();
}

void Data_Mapper::temporal_stuff_plotter(Diff_Eq_Solver *DES) {
        
        this->create_directory();
        this->UI->cout_str("DM-> saving temporal datas in " + this->data_dir);

        this->file.open(this->data_dir + "/temporal_stuff_" + to_string(DES->ite_count + 1) +"_ite.data", ios::out);
        this->write_args();

        this->file << "#" << setw(WIDTH) << "time" 
             <<  setw(WIDTH) << "time_step"
             << setw(WIDTH) << "thrust"
             << setw(WIDTH) << "max pressure"
             << setw(WIDTH) << "max temperature"
             << setw(WIDTH) << "max vol_mass"
             << setw(WIDTH) << "max speed0"
             << setw(WIDTH) << "max speed1"
             << setw(WIDTH) << "max turb_en"
             << setw(WIDTH) << "max turb_dis"
             << endl << endl;

        double time = 0;
        for (int k=0; k<this->arglist_pt->iter_number_solver;k++) {
                this->file << setw(WIDTH) << setprecision(PRECISION) << time
                    << setw(WIDTH) << setprecision(PRECISION) << DES->time_steps[k]
                    << setw(WIDTH) << setprecision(PRECISION) << DES->thrust[k] 
                    << setw(WIDTH) << setprecision(PRECISION) << DES->variables_max.pressure[k]
                    << setw(WIDTH) << setprecision(PRECISION) << DES->variables_max.temperature[k]
                    << setw(WIDTH) << setprecision(PRECISION) << DES->variables_max.vol_mass[k]
                    << setw(WIDTH) << setprecision(PRECISION) << DES->variables_max.speed0[k]
                    << setw(WIDTH) << setprecision(PRECISION) << DES->variables_max.speed1[k]
                    << setw(WIDTH) << setprecision(PRECISION) << DES->variables_max.turb_en[k]
                    << setw(WIDTH) << setprecision(PRECISION) << DES->variables_max.turb_dis[k]
                    << endl;
                time += DES->time_steps[k];
        }
        this->file.close();
}

void show()
{

}
