/* This is nozzle profiler implementation linked to nozzle_profiler.h
 *
 *
 *
 *
*/


// ======== includes ======== 
#include "../include/nozzle_profiler.h"
#include "../include/arg_interpreter.h"
#include "../include/usr_interface.h"
#include "../include/data_mapper.h"
#include <cmath>
// ======== implementation ========


// Fonctions générales appelées depuis l'exterieur

Nozzle_Profiler::Nozzle_Profiler(Usr_Interface *UI, Data_Mapper *DM, arglist_struct *arglist_pt)
{
        this->UI = UI;
        this->DM = DM;
        this->arglist_pt = arglist_pt;
        this->DES = new Diff_Eq_Solver (UI, arglist_pt, &(this->mesh_grid_1), &(this->mesh_grid_2));
        this->create_mesh_grids();
}


void Nozzle_Profiler::profile()
{
        switch ( this->arglist_pt->nozzle_fitting_algo )
        {
                case BRUTAL_FORCE: 
                        this->profile_brutal_force();
                        break;

                case LAGRANGE:
                        this->profile_lagrange();
                        break;

                case SEGMENT:
                        this->profile_segment();
                        break;

                default:
                        throw "NP: Bad profilling algorithm";
                        break;
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

void Nozzle_Profiler::set_init_conditions()
{
        
        float &atmo_press = this->arglist_pt->init_cond.atmosphere_pressure;
        float &atmo_temp = this->arglist_pt->init_cond.atmosphere_temp;
        float &atmo_v_mass = this->arglist_pt->init_cond.atmosphere_vol_mass;
        float &atmo_speed = this->arglist_pt->init_cond.atmosphere_speed;
        float &chbr_press = this->arglist_pt->init_cond.chamber_pressure;
        float &chbr_temp = this->arglist_pt->init_cond.chamber_temp;
        float &chbr_v_mass = this->arglist_pt->init_cond.chamber_vol_mass;
        float &chbr_speed = this->arglist_pt->init_cond.chamber_speed;
        float middle_press, middle_temp, middle_v_mass, middle_speed;
        int middle = this->arglist_pt->y_size - this->arglist_pt->chamber_length -this->arglist_pt->nozzle_length;
        
        register int i,j;

        for(j=0; j<middle; j++){ 
                for(i=0; i<this->arglist_pt->x_size; i++) {
                        this->mesh_grid_1[i][j].pressure = (atmo_press * (this->arglist_pt->y_size -1 - j) + chbr_press * j)/(this->arglist_pt->y_size -1);
                        
                        this->mesh_grid_1[i][j].temperature = (atmo_temp * (this->arglist_pt->y_size -1 - j) + chbr_temp * j)/(this->arglist_pt->y_size -1);
                        
                        this->mesh_grid_1[i][j].vol_mass = (atmo_v_mass * (this->arglist_pt->y_size -1 - j) + chbr_v_mass * j)/(this->arglist_pt->y_size -1);
                        
                        this->mesh_grid_1[i][j].speed[1] = (atmo_speed * (this->arglist_pt->y_size -1 - j) + chbr_speed * j)/(this->arglist_pt->y_size -1);
                }
        }

        middle_press = (atmo_press * (this->arglist_pt->y_size -1 - j) + chbr_press * j)/(this->arglist_pt->y_size -1);
        middle_temp = (atmo_temp * (this->arglist_pt->y_size -1 - j) + chbr_temp * j)/(this->arglist_pt->y_size -1);
        middle_v_mass = (atmo_v_mass * (this->arglist_pt->y_size -1 - j) + chbr_v_mass * j)/(this->arglist_pt->y_size -1);
        middle_speed = (atmo_speed * (this->arglist_pt->y_size -1 - j) + chbr_speed * j)/(this->arglist_pt->y_size -1);

        for (j=middle; j<this->arglist_pt->y_size;j++){
                for(i=0;(not (this->mesh_grid_1[i+1][j].is_wall)); i++){

                        this->mesh_grid_1[i][j].pressure = (middle_press * (this->arglist_pt->y_size -1 - j) 
                                                            + atmo_press * (j-middle))/(this->arglist_pt->y_size -1 - middle);
                        
                        this->mesh_grid_1[i][j].temperature = (middle_temp * (this->arglist_pt->y_size -1 - j) 
                                                            + atmo_temp * (j-middle))/(this->arglist_pt->y_size -1 - middle);
                        
                        this->mesh_grid_1[i][j].vol_mass = (middle_v_mass * (this->arglist_pt->y_size -1 - j) 
                                                            + atmo_v_mass * (j-middle))/(this->arglist_pt->y_size -1 - middle);
                        
                        this->mesh_grid_1[i][j].speed[1] = (middle_speed * (this->arglist_pt->y_size -1 - j) 
                                                            + atmo_speed * (j-middle))/(this->arglist_pt->y_size -1 - middle);
                }

                for(i=i; i<this->arglist_pt->x_size;i++){
                        this->mesh_grid_1[i][j].pressure = (atmo_press * (this->arglist_pt->y_size -1 - j) + chbr_press * j)/(this->arglist_pt->y_size -1);
                        
                        this->mesh_grid_1[i][j].temperature = (atmo_temp * (this->arglist_pt->y_size -1 - j) + chbr_temp * j)/(this->arglist_pt->y_size -1);
                        
                        this->mesh_grid_1[i][j].vol_mass = (atmo_v_mass * (this->arglist_pt->y_size -1 - j) + chbr_v_mass * j)/(this->arglist_pt->y_size -1);
                        
                        this->mesh_grid_1[i][j].speed[1] = (atmo_speed * (this->arglist_pt->y_size -1 - j) + chbr_speed * j)/(this->arglist_pt->y_size -1);

                }
        }
        
        for (i=0; i<this->arglist_pt->x_size; i++) {
                for (j=0; j<this->arglist_pt->y_size; j++) {
                        if (this->mesh_grid_1[i][j].is_wall) {
                                this->mesh_grid_1[i][j].speed[1]=0 ;
                                this->mesh_grid_1[i][j].speed[0]=0 ;
                        }
                }
        }
        
        for (i=0; i<this->arglist_pt->x_size; i++){
                for (j=0; j<this->arglist_pt->y_size; j++){
                        this->mesh_grid_2[i][j].pressure = this->mesh_grid_1[i][j].pressure ;
                        this->mesh_grid_2[i][j].temperature = this->mesh_grid_1[i][j].temperature ;
                        this->mesh_grid_2[i][j].vol_mass = this->mesh_grid_1[i][j].vol_mass ;
                        this->mesh_grid_2[i][j].speed[1] = this->mesh_grid_1[i][j].speed[1] ;

                }
        }
}

void Nozzle_Profiler::save_mesh_grid()
{
        this->DM->create_datafile_from_mesh_grid(&(this->mesh_grid_1));
}


bool Nozzle_Profiler::is_in_x_range(int i)
{
        return (i <= this->arglist_pt->x_size and i >= 0) ;
}
// -- Algo force brute 

void Nozzle_Profiler::init_profile_brutal_force()
{
        int i;       
}

void Nozzle_Profiler::profile_brutal_force()
{
        int i;       
}

void Nozzle_Profiler::one_iteration_brutal_force()
{
        int i;       
}

// -- Algo Lagrange

float Nozzle_Profiler::lagrange_chamber(float X)
{
        float Y=0;
        float Z=1;
        register int i,j;
        for (i=0; i<this->arglist_pt->nozzle_fitting_init_arg.lagrange.nb_pts_chamber; i++){
                for(j=0; j<i; j++){
                      Z *= (X - this->arglist_pt->nozzle_fitting_init_arg.lagrange.abscisses_chamber[j])/
                           (this->arglist_pt->nozzle_fitting_init_arg.lagrange.abscisses_chamber[i] - 
                            this->arglist_pt->nozzle_fitting_init_arg.lagrange.abscisses_chamber[j]);  
                }
                for(j=i+1; j<this->arglist_pt->nozzle_fitting_init_arg.lagrange.nb_pts_chamber; j++){
                      Z *= (X - this->arglist_pt->nozzle_fitting_init_arg.lagrange.abscisses_chamber[j])/
                           (this->arglist_pt->nozzle_fitting_init_arg.lagrange.abscisses_chamber[i] - 
                            this->arglist_pt->nozzle_fitting_init_arg.lagrange.abscisses_chamber[j]);  
                }
                Y += Z;
        }
        return Y;
}

float Nozzle_Profiler::lagrange_nozzle(float X)
{
        float Y=0;
        float Z=1;
        register int i,j;
        for (i=0; i<this->arglist_pt->nozzle_fitting_init_arg.lagrange.nb_pts_nozzle; i++){
                for(j=0; j<i; j++){
                      Z *= (X - this->arglist_pt->nozzle_fitting_init_arg.lagrange.abscisses_nozzle[j])/
                           (this->arglist_pt->nozzle_fitting_init_arg.lagrange.abscisses_nozzle[i] - 
                            this->arglist_pt->nozzle_fitting_init_arg.lagrange.abscisses_nozzle[j]);  
                }
                for(j=i+1; j<this->arglist_pt->nozzle_fitting_init_arg.lagrange.nb_pts_nozzle; j++){
                      Z *= (X - this->arglist_pt->nozzle_fitting_init_arg.lagrange.abscisses_nozzle[j])/
                           (this->arglist_pt->nozzle_fitting_init_arg.lagrange.abscisses_nozzle[i] - 
                            this->arglist_pt->nozzle_fitting_init_arg.lagrange.abscisses_nozzle[j]);  
                }
                Y += Z;
        }
        return Y;
}

void Nozzle_Profiler::init_profile_lagrange()
{
        this->UI->cout_str("NP-> initializing lagrange profile...");
        register int i,j;
        for (j=this->arglist_pt->y_size - this->arglist_pt->chamber_length - this->arglist_pt->nozzle_length; j<this->arglist_pt->y_size - this->arglist_pt->chamber_length; j++){
                
                i = floor(this->lagrange_nozzle(i*this->arglist_pt->space_step));
                if (this->is_in_x_range(i)){
                        this->mesh_grid_1[i][j].is_wall = true;
                        this->mesh_grid_2[i][j].is_wall = true;
                }
                else {throw "NP: impossible to draw initial profile: some points are outside the matrix";}
        }

        for (j=this->arglist_pt->y_size - this->arglist_pt->chamber_length; j<this->arglist_pt->y_size; j++){
                i = floor(this->lagrange_chamber(i*this->arglist_pt->space_step));
                if (this->is_in_x_range(i)) {
                        this->mesh_grid_1[i][j].is_wall = true;
                        this->mesh_grid_2[i][j].is_wall = true;
                }
                else {throw "NP: impossible to draw initial profile: some points are outside the matrix";}
        }
}

void Nozzle_Profiler::profile_lagrange()
{
        this->init_profile_lagrange();
        
        register int i;
        for (i=0; i<this->arglist_pt->iter_number_profiler; i++){ 
                this->one_iteration_lagrange();
        }
}

void Nozzle_Profiler::one_iteration_lagrange()
{
        this->DES->solve();               
}



// fonctions segment


int Nozzle_Profiler::segment(int a_abs, int a_ord, int b_abs, int b_ord, int x_abs)
{
        float x_ord = (float)(a_ord*(b_abs-x_abs) + b_ord*(x_abs-a_abs))/(b_abs-a_abs);
        return floor(x_ord);
}


void Nozzle_Profiler::init_profile_segment()
{
        this->UI->cout_str("NP-> initializing segment profile...");
        register int a,i,j;
        for (a=0; a<this->arglist_pt->nozzle_fitting_init_arg.segment.nb_pts -1 ; a++){
                for (j=this->arglist_pt->nozzle_fitting_init_arg.segment.abscisses[a]; j<this->arglist_pt->nozzle_fitting_init_arg.segment.abscisses[a+1];j++){
                        i = this->arglist_pt->x_size - this->segment(this->arglist_pt->nozzle_fitting_init_arg.segment.abscisses[a],
                                          this->arglist_pt->nozzle_fitting_init_arg.segment.ordinates[a],
                                          this->arglist_pt->nozzle_fitting_init_arg.segment.abscisses[a+1],
                                          this->arglist_pt->nozzle_fitting_init_arg.segment.ordinates[a+1],j);
                        this->mesh_grid_1[i][j].is_wall = true;
                        this->mesh_grid_2[i][j].is_wall = true;
                }
        }

        j = this->arglist_pt->nozzle_fitting_init_arg.segment.abscisses[this->arglist_pt->nozzle_fitting_init_arg.segment.nb_pts -1];
        i = this->arglist_pt->x_size - this->arglist_pt->nozzle_fitting_init_arg.segment.ordinates[this->arglist_pt->nozzle_fitting_init_arg.segment.nb_pts -1];
        
        this->mesh_grid_1[i][j].is_wall = true;
        this->mesh_grid_2[i][j].is_wall = true;
}

void Nozzle_Profiler::profile_segment()
{
        this->init_profile_segment();
        register int i;
        for (i=0; i<this->arglist_pt->iter_number_profiler; i++){
                this->UI->cout_str_no_endl("NP-> iteration nb: "); this->UI->cout_int(i);
                this->one_iteration_segment();
        }
}

void Nozzle_Profiler::one_iteration_segment()
{
        this->set_init_conditions();
        this->DES->solve();
}
