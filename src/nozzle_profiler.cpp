/* This is nozzle profiler implementation linked to nozzle_profiler.h
 *
 *
 *
 *
*/


// ======== includes ======== 
#include "nozzle_profiler.h"
#include "arg_interpreter.h"
#include "usr_interface.h"
#include "data_mapper.h"
#include <gsl/gsl_multimin.h>
#include <cmath>

// ======== implementation ========

#define R 8.314

// Fonctions générales appelées depuis l'exterieur

Nozzle_Profiler::Nozzle_Profiler(Usr_Interface *UI, Data_Mapper *DM, arglist_struct *arglist_pt)
{
        this->UI = UI;
        this->DM = DM;
        this->arglist_pt = arglist_pt;
        this->DES = new Diff_Eq_Solver (UI, this->DM, this, arglist_pt, &(this->mesh_grid_1), &(this->mesh_grid_2));
        this->create_mesh_grids();
}

void Nozzle_Profiler::profile()
{
    switch ( this->arglist_pt->init_cond_type )
    {
        case INIT_GRAD:
             try {this->profile_init_grad();}
             catch (const char *err[]) {throw *err;}
             break;
        case EVOL_CHAMBER:
             try {this->profile_evol_chamber();}
             catch (const char *err[]) {throw *err;}
             break;
    }
}

void Nozzle_Profiler::profile_evol_chamber()
{
        this->init_profile_segment();
        //this-init_profile_constant();
        register int i;
        for (i=0; i<this->arglist_pt->iter_number_profiler; i++){
                this->UI->cout_str_no_endl("NP-> iteration nb: "); this->UI->cout_int(i);
                try {this->one_iteration_evol_chamber();}
                catch (const char *err) {throw *err;}
        }
        this->save_mesh_grid();
}

void Nozzle_Profiler::profile_init_grad()
{
        this->init_profile_segment();
        //this-init_profile_constant();
        register int i;
        for (i=0; i<this->arglist_pt->iter_number_profiler; i++){
                this->UI->cout_str_no_endl("NP-> iteration nb: "); this->UI->cout_int(i);
                try {this->one_iteration_init_grad();}
                catch (const char *err) {throw *err;}
        }
        this->save_mesh_grid();
}


// Algos internes


void Nozzle_Profiler::create_mesh_grids()
{
        this->mesh_grid_1.resize(this->arglist_pt->x_size);
        this->mesh_grid_2.resize(this->arglist_pt->x_size);
        for ( register int i=0; i< this->arglist_pt->x_size; i++){
                this->mesh_grid_1[i].resize(this->arglist_pt->y_size);
                this->mesh_grid_2[i].resize(this->arglist_pt->y_size);
        }
}

void Nozzle_Profiler::set_init_grad()
{
        
        float &mol_mass = this->arglist_pt->mol_mass;
        float &atmo_press = this->arglist_pt->init_cond.atmosphere_pressure;
        float &atmo_temp = this->arglist_pt->init_cond.atmosphere_temp;
        float &atmo_speed = this->arglist_pt->init_cond.atmosphere_speed;
        float &atmo_turb_en = this->arglist_pt->init_cond.atmosphere_turb_en;
        float &atmo_turb_dis = this->arglist_pt->init_cond.atmosphere_turb_dis;
        float atmo_v_mass = atmo_press * mol_mass / (R * atmo_temp);

        float &chbr_press = this->arglist_pt->init_cond.chamber_pressure;
        float &chbr_temp = this->arglist_pt->init_cond.chamber_temp;
        float &chbr_speed = this->arglist_pt->init_cond.chamber_speed;
        float &chbr_turb_en = this->arglist_pt->init_cond.chamber_turb_en;
        float &chbr_turb_dis = this->arglist_pt->init_cond.chamber_turb_dis;
        
        int middle = this->arglist_pt->nozzle_fitting_init_arg.abscisses[0];
        
        
        register int i,j;

        for (j=0; j<middle;j++){
                for (i=0; i<this->arglist_pt->x_size; i++){
                        this->mesh_grid_1[i][j].pressure = atmo_press;
                        this->mesh_grid_1[i][j].temperature = atmo_temp;  
                        this->mesh_grid_1[i][j].vol_mass = atmo_v_mass; 
                        this->mesh_grid_1[i][j].speed[1] = atmo_speed;
                        this->mesh_grid_1[i][j].turb_en = atmo_turb_en;
                        this->mesh_grid_1[i][j].turb_dis = atmo_turb_dis;
                }

        }

        for(j=middle; j<this->arglist_pt->y_size; j++){ 
                
                for(i=0; (( this->mesh_grid_1[i][j].is_wall and this->mesh_grid_1[i+1][j].is_wall) or not(this->mesh_grid_1[i][j].is_wall)); i++){
                        this->mesh_grid_1[i][j].pressure = atmo_press ;
                        this->mesh_grid_1[i][j].temperature = atmo_temp;  
                        this->mesh_grid_1[i][j].vol_mass = atmo_v_mass ;
                        this->mesh_grid_1[i][j].speed[1] = atmo_speed ;
                        this->mesh_grid_1[i][j].turb_en = atmo_turb_en;
                        this->mesh_grid_1[i][j].turb_dis = atmo_turb_dis;
                }
                for(i=i; i<this->arglist_pt->x_size; i++) {
                        this->mesh_grid_1[i][j].pressure = (atmo_press * (this->arglist_pt->y_size -1 - j) 
                                                            + chbr_press * (j - middle)) / (this->arglist_pt->y_size - 1 - middle);

                        this->mesh_grid_1[i][j].temperature = (atmo_temp * (this->arglist_pt->y_size -1 - j) 
                                                               + chbr_temp * (j - middle)) / (this->arglist_pt->y_size - 1 - middle);

                        this->mesh_grid_1[i][j].vol_mass = this->mesh_grid_1[i][j].pressure * mol_mass / (R * this->mesh_grid_1[i][j].temperature);

                        this->mesh_grid_1[i][j].speed[1] = (atmo_speed * (this->arglist_pt->y_size -1 - j) 
                                                            + chbr_speed * (j - middle)) / (this->arglist_pt->y_size - 1 - middle);
                                                            
                        this->mesh_grid_1[i][j].turb_en = (atmo_turb_en * (this->arglist_pt->y_size -1 - j)
                                                           + chbr_turb_en * (j - middle)) / (this->arglist_pt->y_size - 1 - middle);
                        
                        this->mesh_grid_1[i][j].turb_dis = (atmo_turb_dis * (this->arglist_pt->y_size -1 - j)
                                                           + chbr_turb_dis * (j - middle)) / (this->arglist_pt->y_size - 1 - middle);
                }
        }
        
        for (i=0; i<this->arglist_pt->x_size; i++) {
                for (j=0; j<this->arglist_pt->y_size; j++) {
                        if (this->mesh_grid_1[i][j].is_wall) {
                                this->mesh_grid_1[i][j].speed[1]=0 ;
                                this->mesh_grid_1[i][j].speed[0]=0 ;
                                this->mesh_grid_1[i][j].turb_en=0;
                                this->mesh_grid_1[i][j].turb_dis=0;
                        }
                        this->mesh_grid_2[i][j].speed[0] = this->mesh_grid_1[i][j].speed[0];
                        this->mesh_grid_2[i][j].speed[1] = this->mesh_grid_1[i][j].speed[1];
                        this->mesh_grid_2[i][j].temperature = this->mesh_grid_1[i][j].temperature;
                        this->mesh_grid_2[i][j].pressure = this->mesh_grid_1[i][j].pressure;
                        this->mesh_grid_2[i][j].vol_mass = this->mesh_grid_1[i][j].vol_mass;
                        this->mesh_grid_2[i][j].turb_en = this->mesh_grid_1[i][j].turb_en;
                        this->mesh_grid_2[i][j].turb_dis = this->mesh_grid_1[i][j].turb_dis;
                }
        }
}


void Nozzle_Profiler::set_init_cond_evol_chamber()
{
    register int i,j;
    float &mol_mass = this->arglist_pt->mol_mass;
    float &atmo_press = this->arglist_pt->init_cond.atmosphere_pressure;
    float &atmo_temp = this->arglist_pt->init_cond.atmosphere_temp;
    float &atmo_speed = this->arglist_pt->init_cond.atmosphere_speed;
    float &atmo_turb_en = this->arglist_pt->init_cond.atmosphere_turb_en;
    float &atmo_turb_dis = this->arglist_pt->init_cond.atmosphere_turb_dis;
    float atmo_v_mass = atmo_press * mol_mass / (R * atmo_temp);

    for (j=0; j<this->arglist_pt->y_size;j++){
            for (i=0; i<this->arglist_pt->x_size; i++){
                    this->mesh_grid_1[i][j].pressure = atmo_press;
                    this->mesh_grid_1[i][j].temperature = atmo_temp;  
                    this->mesh_grid_1[i][j].vol_mass = atmo_v_mass; 
                    this->mesh_grid_1[i][j].speed[1] = atmo_speed;
                    this->mesh_grid_1[i][j].turb_en = atmo_turb_en;
                    this->mesh_grid_1[i][j].turb_dis = atmo_turb_dis;
            }

    }

    for (i=0; i<this->arglist_pt->x_size; i++) {
            for (j=0; j<this->arglist_pt->y_size; j++) {
                    if (this->mesh_grid_1[i][j].is_wall) {
                            this->mesh_grid_1[i][j].speed[1]=0 ;
                            this->mesh_grid_1[i][j].speed[0]=0 ;
                            this->mesh_grid_1[i][j].turb_en=0;
                            this->mesh_grid_1[i][j].turb_dis=0;
                    }
                    this->mesh_grid_2[i][j].speed[0] = this->mesh_grid_1[i][j].speed[0];
                    this->mesh_grid_2[i][j].speed[1] = this->mesh_grid_1[i][j].speed[1];
                    this->mesh_grid_2[i][j].temperature = this->mesh_grid_1[i][j].temperature;
                    this->mesh_grid_2[i][j].pressure = this->mesh_grid_1[i][j].pressure;
                    this->mesh_grid_2[i][j].vol_mass = this->mesh_grid_1[i][j].vol_mass;
                    this->mesh_grid_2[i][j].turb_en = this->mesh_grid_1[i][j].turb_en;
                    this->mesh_grid_2[i][j].turb_dis = this->mesh_grid_1[i][j].turb_dis;
            }
    }
}


void Nozzle_Profiler::save_mesh_grid()
{
        this->DM->create_datafile_from_mesh_grid(this->DES);
}

void Nozzle_Profiler::set_wall(int i,int j)
{
    if (i > 0 and j < this->arglist_pt->y_size){ this->mesh_grid_1[i][j].is_wall = true; this->mesh_grid_2[i][j].is_wall = true;}
    if (i-1 > 0 and j < this->arglist_pt->y_size){ this->mesh_grid_1[i-1][j].is_wall = true; this->mesh_grid_2[i-1][j].is_wall = true;}
    if (i > 0 and j+1 < this->arglist_pt->y_size){ this->mesh_grid_1[i][j+1].is_wall = true; this->mesh_grid_2[i][j+1].is_wall = true;}
    if (i-1 > 0 and j+1 < this->arglist_pt->y_size){ this->mesh_grid_1[i-1][j+1].is_wall = true; this->mesh_grid_2[i-1][j+1].is_wall = true;}

}

bool Nozzle_Profiler::is_in_x_range(int i)
{
        return (i <= this->arglist_pt->x_size and i >= 0) ;
}


int segment(int a_abs, int a_ord, int b_abs, int b_ord, int x_abs)
{
        float x_ord = (float)(a_ord*(b_abs-x_abs) + b_ord*(x_abs-a_abs))/(b_abs-a_abs);
        return floor(x_ord);
}


void Nozzle_Profiler::init_profile_segment()
{
        this->UI->cout_str("NP-> initializing segment profile...");
        register int a,j,i_new, i_mem, i;
        i_mem = this->arglist_pt->x_size - this->arglist_pt->nozzle_fitting_init_arg.ordinates[0];

        for (a=0; a<this->arglist_pt->nozzle_fitting_init_arg.nb_pts -1; a++){
                for (j=this->arglist_pt->nozzle_fitting_init_arg.abscisses[a]; j<this->arglist_pt->nozzle_fitting_init_arg.abscisses[a+1];j++){
                        i_new = this->arglist_pt->x_size - segment(this->arglist_pt->nozzle_fitting_init_arg.abscisses[a],
                                          this->arglist_pt->nozzle_fitting_init_arg.ordinates[a],
                                          this->arglist_pt->nozzle_fitting_init_arg.abscisses[a+1],
                                          this->arglist_pt->nozzle_fitting_init_arg.ordinates[a+1],j);
                        this->set_wall(i_new, j);
                        
                        if (i_new > i_mem){
                                for (i=i_new; i>i_mem; i--){
                                        this->set_wall(i,j-1);
                                }
                        }
                        if (i_new < i_mem){
                                for (i=i_new; i<i_mem; i++){
                                        this->set_wall(i,j-1);
                                }
                        }

                        i_mem = i_new;
                }
        }

        j = this->arglist_pt->nozzle_fitting_init_arg.abscisses[this->arglist_pt->nozzle_fitting_init_arg.nb_pts -1];
        i_new = this->arglist_pt->x_size - this->arglist_pt->nozzle_fitting_init_arg.ordinates[this->arglist_pt->nozzle_fitting_init_arg.nb_pts -1];
        
        this->set_wall(i_new,j);
        this->set_wall(i_new,j);

        if (i_new > i_mem){
                for (i=i_new; i>i_mem; i--){
                        this->set_wall(i,j-1);
                }
        }
        if (i_new < i_mem){
                for (i=i_new; i<i_mem; i++){
                        this->set_wall(i,j-1);
                }
        }


}


void Nozzle_Profiler::init_profile_constant()
{
        register int i;
        for (i=this->arglist_pt->y_size/2; i<this->arglist_pt->y_size; i++) {
            this->set_wall(this->arglist_pt->x_size/2,i);
        }
}


void Nozzle_Profiler::one_iteration_init_grad()
{
        this->set_init_grad();
        try {this->DES->solve();}
        catch (const char *err[]) {throw *err;}
}


void Nozzle_Profiler::one_iteration_evol_chamber()
{
        this->set_init_cond_evol_chamber();
        try {this->DES->solve();}
        catch (const char *err[]) {throw *err;}
}

void Nozzle_Profiler::update_chamber_cond()
{
        float &mol_mass = this->arglist_pt->mol_mass;
        float &atmo_pressure = this->arglist_pt->init_cond.atmosphere_pressure;
        float &atmo_temperature = this->arglist_pt->init_cond.atmosphere_temp;
        float atmo_vol_mass = atmo_pressure * mol_mass / (R * atmo_temperature);
        float &atmo_speed = this->arglist_pt->init_cond.atmosphere_speed;
        float &atmo_turb_en = this->arglist_pt->init_cond.atmosphere_turb_en;
        float &atmo_turb_dis = this->arglist_pt->init_cond.atmosphere_turb_dis;

        float &chbr_pressure = this->arglist_pt->init_cond.chamber_pressure;
        float &chbr_temperature = this->arglist_pt->init_cond.chamber_temp;
        float chbr_vol_mass = chbr_pressure * mol_mass / (R * chbr_temperature);
        float &chbr_speed = this->arglist_pt->init_cond.chamber_speed;
        float &chbr_turb_en = this->arglist_pt->init_cond.chamber_turb_en;
        float &chbr_turb_dis = this->arglist_pt->init_cond.chamber_turb_dis;

        register int i;
        for (i=this->arglist_pt->x_size - 1; not (this->mesh_grid_2[i][this->arglist_pt->y_size -1 ].is_wall); i--){
                this->mesh_grid_2[i][this->arglist_pt->y_size -1 ].pressure = 
                    (chbr_pressure * this->DES->ite_count + atmo_pressure * (this->arglist_pt->iter_number_evol_chamber - this->DES->ite_count) ) 
                    / this->arglist_pt->iter_number_evol_chamber;
                this->mesh_grid_2[i][this->arglist_pt->y_size -1 ].temperature = 
                    (chbr_temperature * this->DES->ite_count + atmo_temperature * (this->arglist_pt->iter_number_evol_chamber - this->DES->ite_count) ) 
                    / this->arglist_pt->iter_number_evol_chamber;
                this->mesh_grid_2[i][this->arglist_pt->y_size -1 ].vol_mass = 
                    (chbr_vol_mass * this->DES->ite_count + atmo_vol_mass * (this->arglist_pt->iter_number_evol_chamber - this->DES->ite_count) ) 
                    / this->arglist_pt->iter_number_evol_chamber;
                this->mesh_grid_2[i][this->arglist_pt->y_size -1 ].speed[1] = 
                    (chbr_speed * this->DES->ite_count + atmo_speed * (this->arglist_pt->iter_number_evol_chamber - this->DES->ite_count) ) 
                    / this->arglist_pt->iter_number_evol_chamber;
                this->mesh_grid_2[i][this->arglist_pt->y_size -1 ].turb_en = 
                    (chbr_turb_en * this->DES->ite_count + atmo_turb_en * (this->arglist_pt->iter_number_evol_chamber - this->DES->ite_count) ) 
                    / this->arglist_pt->iter_number_evol_chamber;
                this->mesh_grid_2[i][this->arglist_pt->y_size -1 ].turb_dis = 
                    (chbr_turb_dis * this->DES->ite_count + atmo_turb_dis * (this->arglist_pt->iter_number_evol_chamber - this->DES->ite_count) ) 
                    / this->arglist_pt->iter_number_evol_chamber;
        }
}
