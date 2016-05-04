/* This is differential equation solver implementation linked to diff_eq_solver.h
 *
 *
 *
 *
*/


// ======== includes ======== 
#include "diff_eq_solver.h"
#include "arg_interpreter.h"
#include "nozzle_profiler.h"
#include "data_mapper.h"
#include "usr_interface.h"
#include <cmath>
#include <iostream>

#define R 8.314 //constante des gaz parfaits
//toutes les constantes relatives au modèle k-epsilon
#define sig_k 1
#define sig_e 1.3
#define c_e_1 1.45
#define c_e_2 1.92
#define c_mu 0.09
#define sig_t 1.3

#define mesh_grid_1 (*(this->mesh_grid_pt1))
#define mesh_grid_2 (*(this->mesh_grid_pt2))
#define time_step this->time_steps[this->time_steps.back()]

// ======== implementation ========

//le contruscteur : on crée un vecteur thrust qui va contenir la liste des valeurs de la poussée au cours du temps
Diff_Eq_Solver::Diff_Eq_Solver(Usr_Interface *UI,Data_Mapper *DM, Nozzle_Profiler * NP, arglist_struct *arglist_pt, mesh_grid_t *mesh_grid_pt1, mesh_grid_t *mesh_grid_pt2) {
        this->UI = UI;
        this->DM = DM;
        this->NP = NP;
        this->arglist_pt = arglist_pt;
        this->mesh_grid_pt1 = mesh_grid_pt1;
        this->mesh_grid_pt2 = mesh_grid_pt2;
        this->thrust.resize( this->arglist_pt->iter_number_solver );
        this->time_steps.resize( this->arglist_pt->iter_number_solver );
        this->threads.resize(this->arglist_pt->nb_of_threads);
      
        this->variables_max.pressure_loc.resize(this->arglist_pt->nb_of_threads);
        this->variables_max.temperature_loc.resize(this->arglist_pt->nb_of_threads);
        this->variables_max.vol_mass_loc.resize(this->arglist_pt->nb_of_threads);
        this->variables_max.speed0_loc.resize(this->arglist_pt->nb_of_threads);
        this->variables_max.speed1_loc.resize(this->arglist_pt->nb_of_threads);
        this->variables_max.turb_en_loc.resize(this->arglist_pt->nb_of_threads);
        this->variables_max.turb_dis_loc.resize(this->arglist_pt->nb_of_threads);

        this->variables_max.pressure.resize(this->arglist_pt->iter_number_solver);
        this->variables_max.temperature.resize(this->arglist_pt->iter_number_solver);
        this->variables_max.vol_mass.resize(this->arglist_pt->iter_number_solver);
        this->variables_max.speed0.resize(this->arglist_pt->iter_number_solver);
        this->variables_max.speed1.resize(this->arglist_pt->iter_number_solver);
        this->variables_max.turb_en.resize(this->arglist_pt->iter_number_solver);
        this->variables_max.turb_dis.resize(this->arglist_pt->iter_number_solver);
}


//Fonction qui cherche le maximum de chaque valeur sur le mesh grid

void Diff_Eq_Solver::find_all_max()
{
    register int i, j, t, i_min, i_max, slice;
    i_min = 1;
    slice = (this->arglist_pt->x_size - 1)/ this->arglist_pt->nb_of_threads; // taille d'une tranche horizontale
    i_max = i_min + slice;

    for (t = 0 ; t < this->arglist_pt->nb_of_threads - 1; t++){ 
        this->threads[t] = thread (&Diff_Eq_Solver::partial_find_all_max, this, i_min, i_max, t); // on envoie sur les coeurs
        i_min = i_max;
        i_max += slice;
    }

    this->threads[this->arglist_pt->nb_of_threads - 1] = 
        thread (&Diff_Eq_Solver::partial_find_all_max, this, i_min, this->arglist_pt->x_size-1, this->arglist_pt->nb_of_threads-1); // le dernier se tape le rab de la division euclidienne

    for (auto &thrd : this->threads) thrd.join();
}


void Diff_Eq_Solver::partial_find_all_max(int i_min, int i_max, int t)
{
    register int i,j;
    register data_t max_p, max_t, max_v, max_s0, max_s1, max_te, max_td
        , val_p, val_t, val_v, val_s0, val_s1, val_te, val_td;

    max_p = this->get_pressure(i_min, 0);
    max_t = this->get_temperature(i_min, 0);
    max_v = this->get_vol_mass(i_min,0);
    max_s0 = this->get_speed0(i_min, 0);
    max_s1 = this->get_speed1(i_min , 0);
    max_te = this->get_turb_en(i_min, 0);
    max_td = this->get_turb_dis(i_min, 0);
    

    for (i=i_min; i<i_max; i++){
        for (j=0; j<this->arglist_pt->y_size; j++){
            val_p = this->get_pressure(i,j);
            if (val_p > max_p) max_p = val_p;

            val_t = this->get_temperature(i,j);
            if (val_t > max_t) max_t = val_t;

            val_v = this->get_vol_mass(i,j);
            if (val_v > max_v) max_v = val_v;
            
            val_s0 = this->get_speed0(i,j);
            if (val_s0 > max_s0) max_s0 = val_s0;
            
            val_s1 = this->get_speed1(i,j);
            if (val_s1 > max_s1) max_s1 = val_s1;
            
            val_te = this->get_turb_en(i,j);
            if (val_te > max_te) max_te = val_te;

            val_td = this->get_turb_dis(i,j);
            if (val_td > max_td) max_td = val_td;
        }
    }

    this->variables_max.pressure_loc[t] = max_p;
    this->variables_max.temperature_loc[t] = max_t;
    this->variables_max.vol_mass_loc[t] = max_v;
    this->variables_max.speed0_loc[t] = max_s0;
    this->variables_max.speed1_loc[t] = max_s1;
    this->variables_max.turb_en_loc[t] = max_te;
    this->variables_max.turb_dis_loc[t] = max_td;
}

// avoir la valeur absolu d'une variable d'une case

bool is_normal(data_t d) {return d == d;}

data_t Diff_Eq_Solver::get_pressure(int i, int j)
{
    if ( is_normal(mesh_grid_2[i][j].pressure)){ 
        return abs(mesh_grid_2[i][j].pressure) ;
    }
    else{ 
        this->UI->space();
//        this->UI->cout_float(mesh_grid_2[i][j].pressure);
        this->mtx.lock();
        this->DM->create_datafile_from_mesh_grid(this);
        this->DM->temporal_stuff_plotter(this); 
        throw "DES-> nan values appeared (pressure)";
        this->mtx.unlock();
    }
}

data_t Diff_Eq_Solver::get_temperature(int i, int j)
{
    if ( is_normal(mesh_grid_2[i][j].temperature) ) { 
        return abs(mesh_grid_2[i][j].temperature) ;
    }
    else{        
        this->UI->space();
//        this->UI->cout_float(mesh_grid_2[i][j].temperature);
        this->mtx.lock();
        this->DM->create_datafile_from_mesh_grid(this);
        this->DM->temporal_stuff_plotter(this); 
        throw "DES-> nan values appeared (temperature)";
        this->mtx.unlock();
    }
}

data_t Diff_Eq_Solver::get_vol_mass(int i, int j)
{
    if ( is_normal(mesh_grid_2[i][j].vol_mass) ) {
        return abs(mesh_grid_2[i][j].vol_mass) ;
    }
    else{
        this->UI->space();
//        this->UI->cout_float(mesh_grid_2[i][j].vol_mass);
        this->mtx.lock();
        this->DM->create_datafile_from_mesh_grid(this); 
        this->DM->temporal_stuff_plotter(this); 
        throw "DES-> nan values appeared (vol_mass)";
        this->mtx.unlock();
    }
}

data_t Diff_Eq_Solver::get_speed0(int i, int j)
{
    if ( is_normal(mesh_grid_2[i][j].speed[0]) ) { 
        return abs(mesh_grid_2[i][j].speed[0]) ;
    }
    else{
        this->UI->space();
//        this->UI->cout_float(mesh_grid_2[i][j].speed[0]);
        this->mtx.lock();
        this->DM->create_datafile_from_mesh_grid(this);
        this->DM->temporal_stuff_plotter(this); 
        throw "DES-> nan values appeared (speed0)";
        this->mtx.unlock();
    }
}

data_t Diff_Eq_Solver::get_speed1(int i, int j)
{
    if ( is_normal(mesh_grid_2[i][j].speed[1]) ) { 
        return abs(mesh_grid_2[i][j].speed[1]) ;
    }
    else{
        this->UI->space();
//        this->UI->cout_float(mesh_grid_2[i][j].speed[1]);
        this->mtx.lock();
        this->DM->create_datafile_from_mesh_grid(this);
        this->DM->temporal_stuff_plotter(this); 
        throw "DES-> nan values appeared (speed1)";
        this->mtx.unlock();
    }
}

data_t Diff_Eq_Solver::get_turb_en(int i, int j)
{
    if ( is_normal(mesh_grid_2[i][j].turb_en) ) { 
        return abs(mesh_grid_2[i][j].turb_en) ;
    }
    else{
        this->UI->space();
//        this->UI->cout_float(mesh_grid_2[i][j].turb_en);
        this->mtx.lock();
        this->DM->create_datafile_from_mesh_grid(this);
        this->DM->temporal_stuff_plotter(this); 
        throw "DES-> nan values appeared (turb_en)";
        this->mtx.unlock();
    }
}

data_t Diff_Eq_Solver::get_turb_dis(int i, int j)
{
    if ( is_normal(mesh_grid_2[i][j].turb_dis)) { 
        return abs(mesh_grid_2[i][j].turb_dis);
    }
    else {
        this->UI->space();
//        this->UI->cout_float(mesh_grid_2[i][j].turb_dis);
        this->mtx.lock();
        this->DM->create_datafile_from_mesh_grid(this); 
        this->DM->temporal_stuff_plotter(this); 
        throw "DES-> nan values appeared (turb_dis)";
        this->mtx.unlock();
    }
}

double vec_max(vector<data_t> * vec)
{
    double max ,val;
    max = (*vec)[0];
    register int i;
    for (i = 0; i<vec->size(); i++){
        val = (*vec)[i];
        if ( val > max ) max = val;
    }
    return max;
}

void Diff_Eq_Solver::calc_all_max(int i)
{
    this->variables_max.pressure[i] = vec_max(&this->variables_max.pressure_loc);
    this->variables_max.temperature[i] = vec_max(&this->variables_max.temperature_loc);
    this->variables_max.vol_mass[i] = vec_max(&this->variables_max.vol_mass_loc);
    this->variables_max.speed0[i] = vec_max(&this->variables_max.speed0_loc);
    this->variables_max.speed1[i] = vec_max(&this->variables_max.speed1_loc);
    this->variables_max.turb_en[i] = vec_max(&this->variables_max.turb_en_loc);
    this->variables_max.turb_dis[i] = vec_max(&this->variables_max.turb_dis_loc);
}


//fonction qui calcule le carré de la vitesse, le troisième argument vaut 1:mesh_grid_1 2:mesh_grid_2
data_t Diff_Eq_Solver::speed2(int i, int j, mesh_grid_t * mesh_grid_pt) {
        return(pow((*mesh_grid_pt)[i][j].speed[0],2) + pow((*mesh_grid_pt)[i][j].speed[1],2));
}

//la pression totale est rho*e+P
data_t Diff_Eq_Solver::pres_tot_PG(int i, int j) {
        if (mesh_grid_1[i][j].is_wall) {
        	return(0);
        }
        else {
        	return(
        	mesh_grid_1[i][j].vol_mass*en_tot_PG(i,j)+mesh_grid_1[i][j].pressure
        	);
        }
}

data_t Diff_Eq_Solver::pres_tot_VDW(int i, int j) {
	if (mesh_grid_1[i][j].is_wall) {
        	return(0);
        }
        else {
        	return(
        	mesh_grid_1[i][j].vol_mass*en_tot_VDW(i,j)+mesh_grid_1[i][j].pressure
        	);
        }
}

//l'énergie totale est u+e_c
data_t Diff_Eq_Solver::en_tot_PG(int i, int j) {
        if (mesh_grid_1[i][j].is_wall) {
        	return(0);
        }
        else {
        	return(
        	5.0/2.0*R/this->arglist_pt->mol_mass*mesh_grid_1[i][j].temperature + 1.0/2.0*speed2(i,j,this->mesh_grid_pt1)
        	);
        }
}

data_t Diff_Eq_Solver::en_tot_VDW(int i, int j) {
        if (mesh_grid_1[i][j].is_wall) {
        	return(0);
        }
        else {
        	return(
        	5.0/2.0*R*mesh_grid_1[i][j].temperature/this->arglist_pt->mol_mass + 1.0/2.0*speed2(i,j,this->mesh_grid_pt1) - this->arglist_pt->VDW_a_coef*mesh_grid_1[i][j].vol_mass/pow(this->arglist_pt->mol_mass,2)
        	);
        }
}

//divergence de rho*v en cartésiennes
data_t Diff_Eq_Solver::diver_rhov_c(int i, int j) {
        return(
             deriv_x_rhovx(i,j)+deriv_y_rhovy(i,j)
        );
}

//dérivée de la température par rapport à x et y
data_t Diff_Eq_Solver::deriv_x_temp(int i, int j) {
        return(
        (1-mesh_grid_1[i+1][j].is_wall)*(1-mesh_grid_1[i-1][j].is_wall)*(mesh_grid_1[i+1][j].temperature-mesh_grid_1[i-1][j].temperature)
        / (2*this->arglist_pt->space_step)
        );
}

data_t Diff_Eq_Solver::deriv_y_temp(int i, int j) {
        return(
        (1-mesh_grid_1[i][j+1].is_wall)*(1-mesh_grid_1[i][j-1].is_wall)*(mesh_grid_1[i][j+1].temperature-mesh_grid_1[i][j-1].temperature)
        / (2*this->arglist_pt->space_step)
        );
}

//dérivée seconde de la température par rapport à x et y. On l'écrit en éléments finis en faisant gaffe que le grad de T au mur est nul (continuité)
data_t Diff_Eq_Solver::deriv2_x_temp(int i, int j) {
        return(
        ((1-mesh_grid_1[i+1][j].is_wall)*(mesh_grid_1[i+1][j].temperature-mesh_grid_1[i][j].temperature)-(1-mesh_grid_1[i-1][j].is_wall)*(mesh_grid_1[i][j].temperature-mesh_grid_1[i-1][j].temperature))
        / (pow(this->arglist_pt->space_step,2))
        );
}

data_t Diff_Eq_Solver::deriv2_y_temp(int i, int j) {
        return(
        ((1-mesh_grid_1[i][j+1].is_wall)*(mesh_grid_1[i][j+1].temperature-mesh_grid_1[i][j].temperature)-(1-mesh_grid_1[i][j-1].is_wall)*(mesh_grid_1[i][j].temperature-mesh_grid_1[i][j-1].temperature))
        / (pow(this->arglist_pt->space_step,2))
        );
}

//dérivée de la pression par rapport à x et y
data_t Diff_Eq_Solver::deriv_x_pres(int i, int j) {
        return(
        (1-mesh_grid_1[i+1][j].is_wall)*(1-mesh_grid_1[i-1][j].is_wall)*(mesh_grid_1[i+1][j].pressure-mesh_grid_1[i-1][j].pressure)
        / (this->arglist_pt->space_step*2)
        );
}

data_t Diff_Eq_Solver::deriv_y_pres(int i, int j) {
        return(
        (1-mesh_grid_1[i][j+1].is_wall)*(1-mesh_grid_1[i][j-1].is_wall)*(mesh_grid_1[i][j+1].pressure-mesh_grid_1[i][j-1].pressure)
        / (this->arglist_pt->space_step*2)
        );
}

//dérivée par rapport à x de rho*v_x
data_t Diff_Eq_Solver::deriv_x_rhovx(int i, int j) {
        return(
        (mesh_grid_1[i+1][j].vol_mass*mesh_grid_1[i+1][j].speed[0]-mesh_grid_1[i-1][j].vol_mass*mesh_grid_1[i-1][j].speed[0])
        / (2*this->arglist_pt->space_step)
        );
}

//dérivée par rapport à y de rho*v_x
data_t Diff_Eq_Solver::deriv_y_rhovx(int i, int j) {
        return(
        (mesh_grid_1[i][j+1].vol_mass*mesh_grid_1[i][j+1].speed[0]-mesh_grid_1[i][j-1].vol_mass*mesh_grid_1[i][j-1].speed[0])
        / (this->arglist_pt->space_step*2)
        );
}

//dérivée par rapport à x de rho*v_y
data_t Diff_Eq_Solver::deriv_x_rhovy(int i, int j) {
        return(
        (mesh_grid_1[i+1][j].vol_mass*mesh_grid_1[i+1][j].speed[1]-mesh_grid_1[i-1][j].vol_mass*mesh_grid_1[i-1][j].speed[1])
        / (this->arglist_pt->space_step*2)
        );
}

//dérivée par rapport à y de rho*v_y
data_t Diff_Eq_Solver::deriv_y_rhovy(int i, int j) {
        return(
        (mesh_grid_1[i][j+1].vol_mass*mesh_grid_1[i][j+1].speed[1]-mesh_grid_1[i][j-1].vol_mass*mesh_grid_1[i][j-1].speed[1])
        / (this->arglist_pt->space_step*2)
        );
}

//dérivée par rapport à x (y) de la pression_totale*v_x (v_y), pour le gaz parfait
data_t Diff_Eq_Solver::deriv_x_prestotPGvx(int i, int j) {
        return(
        (mesh_grid_1[i+1][j].speed[0]*pres_tot_PG(i+1,j)-mesh_grid_1[i-1][j].speed[0]*pres_tot_PG(i-1,j))
        / (this->arglist_pt->space_step*2)
        );
}

data_t Diff_Eq_Solver::deriv_y_prestotPGvy(int i, int j) {
        return(
        (mesh_grid_1[i][j+1].speed[1]*pres_tot_PG(i,j+1)-mesh_grid_1[i][j-1].speed[1]*pres_tot_PG(i,j-1))
        / (this->arglist_pt->space_step*2)
        );
}

//dérivée par rapport à r (coord cyl) de la pression_totale*r*v_r
data_t Diff_Eq_Solver::deriv_r_prestotPGrvr(int i, int j) {
        return(
        -(r(i+1)*mesh_grid_1[i+1][j].speed[0]*pres_tot_PG(i+1,j)-r(i-1)*mesh_grid_1[i-1][j].speed[0]*pres_tot_PG(i-1,j))
        / (this->arglist_pt->space_step*2)
        );
}

//dérivée par rapport à r (coord cyl) de la température
data_t Diff_Eq_Solver::deriv_r_temp(int i, int j) {
        return(
        (1-mesh_grid_1[i+1][j].is_wall)*(1-mesh_grid_1[i-1][j].is_wall)*(mesh_grid_1[i-1][j].temperature-mesh_grid_1[i+1][j].temperature)
        / (this->arglist_pt->space_step*2)
        );
}

//dérivée par rapport à x (y) de la pression_totale*v_x (v_y), pour le gaz de Van der Waals
data_t Diff_Eq_Solver::deriv_x_prestotVDWvx(int i, int j) {
        return(
        (mesh_grid_1[i+1][j].speed[0]*pres_tot_VDW(i+1,j)-mesh_grid_1[i-1][j].speed[0]*pres_tot_VDW(i-1,j))
        / (this->arglist_pt->space_step*2)
        );
}

data_t Diff_Eq_Solver::deriv_y_prestotVDWvy(int i, int j) {
        return(
        (mesh_grid_1[i][j+1].speed[1]*pres_tot_VDW(i,j+1)-mesh_grid_1[i][j-1].speed[1]*pres_tot_VDW(i,j-1))
        / (this->arglist_pt->space_step*2)
        );
}

//dérivée par rapport à r (coord cyl) de la pression_totale*r*v_r, pour le gaz de Van der Waals
data_t Diff_Eq_Solver::deriv_r_prestotVDWrvr(int i, int j) {
        return(
        -(r(i+1)*mesh_grid_1[i+1][j].speed[0]*pres_tot_VDW(i+1,j)-r(i-1)*mesh_grid_1[i-1][j].speed[0]*pres_tot_VDW(i-1,j))
        / (this->arglist_pt->space_step*2)
        );
}

//définition de la coordonnée cylindrique r (=-i en fait)
data_t Diff_Eq_Solver::r(int i) {
        return(this->arglist_pt->space_step*(this->arglist_pt->x_size-i));
}

//dérivée par rapport à r de r*rho*v_r (cyl)
data_t Diff_Eq_Solver::deriv_r_rrhovr(int i, int j) {
        return(
        (r(i-1)*mesh_grid_1[i-1][j].vol_mass*mesh_grid_1[i-1][j].speed[0]-r(i+1)*mesh_grid_1[i+1][j].vol_mass*mesh_grid_1[i+1][j].speed[0])
        / (this->arglist_pt->space_step*2)
        );
}

//toutes les dérivées de la vitesse en cartésiennes
data_t Diff_Eq_Solver::deriv_x_vx(int i, int j) {
	return(
	(mesh_grid_1[i+1][j].speed[0]-mesh_grid_1[i-1][j].speed[0])
	/ (this->arglist_pt->space_step*2)
	);
}

data_t Diff_Eq_Solver::deriv_y_vx(int i, int j) {
	return(
	(mesh_grid_1[i][j+1].speed[0]-mesh_grid_1[i][j-1].speed[0])
	/ (this->arglist_pt->space_step*2)
	);
}

data_t Diff_Eq_Solver::deriv_x_vy(int i, int j) {
	return(
	(mesh_grid_1[i+1][j].speed[1]-mesh_grid_1[i-1][j].speed[1])
	/ (this->arglist_pt->space_step*2)
	);
}

data_t Diff_Eq_Solver::deriv_y_vy(int i, int j) {
	return(
	(mesh_grid_1[i][j+1].speed[1]-mesh_grid_1[i][j-1].speed[1])
	/ (this->arglist_pt->space_step*2)
	);
}

//le "strain" : c'est le déplacement utilisé dans le calcul de la contrainte (ici en cartésiennes)
//attention que le déplacement est symétrique, strain_xy=strain_yx
data_t Diff_Eq_Solver::strain_xy(int i, int j) {
	if (mesh_grid_1[i][j].is_wall) {
        	return(0);
        }
        else {
		return(1.0/2.0*(deriv_y_vx(i,j)+deriv_x_vy(i,j)));
        }
}

data_t Diff_Eq_Solver::strain_xx(int i, int j) {
	if (mesh_grid_1[i][j].is_wall) {
        	return(0);
        }
        else {
		return(deriv_x_vx(i,j));
        }
}

data_t Diff_Eq_Solver::strain_yy(int i, int j) {
	if (mesh_grid_1[i][j].is_wall) {
        	return(0);
        }
	else {
		return(deriv_y_vy(i,j));
	}
}

//les dérivées de la contrainte moléculaire pour le modèle non turbulent
//attention : double dérivée peut sortir du tableau, si ça sort on met 0
data_t Diff_Eq_Solver::deriv_y_tauxy(int i, int j) {
	if (mesh_grid_1[i][j].is_wall) {
        	return(0);
        }
	if (is_in(i,j+2) && is_in(i,j-2)) {
		return(
		(mol_stress_xy(i,j+1)-mol_stress_xy(i,j-1))
		/ (this->arglist_pt->space_step*2)
		);
	}
	else {
		return(0);
	}
}

data_t Diff_Eq_Solver::deriv_x_tauxx(int i, int j) {
	if (mesh_grid_1[i][j].is_wall) {
        	return(0);
        }
	if (is_in(i+2,j) && is_in(i-2,j)) {
		return(
		(mol_stress_xx(i+1,j)-mol_stress_xx(i-1,j))
		/ (this->arglist_pt->space_step*2)
		);
	}
	else {
		return(0);
	}
}

data_t Diff_Eq_Solver::deriv_y_tauyy(int i, int j) {
	if (mesh_grid_1[i][j].is_wall) {
        	return(0);
        }
	if (is_in(i,j+2) && is_in(i,j-2)) {
		return(
		(mol_stress_yy(i,j+1)-mol_stress_yy(i,j-1))
		/(this->arglist_pt->space_step*2)
		);
	}
	else {
		return(0);
	}
}

data_t Diff_Eq_Solver::deriv_x_tauyx(int i, int j) {
	if (mesh_grid_1[i][j].is_wall) {
        	return(0);
        }
	if (is_in(i+2,j) && is_in(i-2,j)) {
		return(
		(mol_stress_xy(i+1,j)-mol_stress_xy(i-1,j))
		/ (this->arglist_pt->space_step*2)
		);
	}
	else {
		return(0);
	}
}

//la divergence de la contrainte*vitesse (en cartésiennes non turbulentes)
data_t Diff_Eq_Solver::diver_vtau(int i, int j) {
	return(
	deriv_x_vtaux(i,j)+deriv_y_vtauy(i,j)
	);
}

//la contrainte*vitesse (en cartésiennes non turbulentes)
data_t Diff_Eq_Solver::vtaux(int i, int j) {
	return(
	mesh_grid_1[i][j].speed[0]*mol_stress_xx(i,j)+mesh_grid_1[i][j].speed[1]*mol_stress_xy(i,j)
	);
}

data_t Diff_Eq_Solver::vtauy(int i, int j) {
	return(
	mesh_grid_1[i][j].speed[0]*mol_stress_xy(i,j)+mesh_grid_1[i][j].speed[1]*mol_stress_yy(i,j)
	);
}

//les dérivées de contrainte*vitesse (en cartésiennes non tubulentes)
//attention à ne pas sortir du tableau
data_t Diff_Eq_Solver::deriv_x_vtaux(int i, int j) {
	if (mesh_grid_1[i][j].is_wall) {
        	return(0);
        }
	if (is_in(i+2,j) && is_in(i-2,j)) {
		return(
		(vtaux(i+1,j)-vtaux(i-1,j))
		/ (this->arglist_pt->space_step*2)
		);
	}
	else {
		return(0);
	}
}

data_t Diff_Eq_Solver::deriv_y_vtauy(int i, int j) {
	if (mesh_grid_1[i][j].is_wall) {
        	return(0);
	}
	if (is_in(i,j+2) && is_in(i,j-2)) {
		return(
		(vtauy(i,j+1)-vtauy(i,j-1))
		/ (this->arglist_pt->space_step*2)
		);
	}
	else {
		return(0);
	}
}

//définition de la viscosité turbulente
data_t Diff_Eq_Solver::mu_t(int i, int j) {
	if (not(mesh_grid_1[i][j].is_wall)) {
	return( c_mu*mesh_grid_1[i][j].vol_mass*pow(mesh_grid_1[i][j].turb_en,2)/mesh_grid_1[i][j].turb_dis*exp(-3.4/pow(1.0+0.02*Ret(i,j),2)));
	}
	else { return(0); }
}

//définition de la conductivité thermique turbulente
data_t Diff_Eq_Solver::lambda_t(int i, int j) {
	return(
	5.0/2.0*R/this->arglist_pt->mol_mass*mu_t(i,j)/sig_t
	);
}

//définition du nombre de Reynolds turbulent
data_t Diff_Eq_Solver::Ret(int i, int j) {
	if (not(mesh_grid_1[i][j].is_wall) && is_in(i,j)) {
	return(
	mesh_grid_1[i][j].vol_mass*pow(mesh_grid_1[i][j].turb_en,2)/this->arglist_pt->dyn_visc/mesh_grid_1[i][j].turb_dis
	);
	}
	else {
		return(0);
	}
}

//définition du flux de chaleur turbulent (en cartésiennes)
data_t Diff_Eq_Solver::heat_flux_x_turb(int i, int j) {
	if (mesh_grid_1[i][j].is_wall) {
        	return(0);
        }
	else {
		return(-(this->arglist_pt->lambda+lambda_t(i,j))*deriv_x_temp(i,j));
	}
}

data_t Diff_Eq_Solver::heat_flux_y_turb(int i, int j) {
	if (mesh_grid_1[i][j].is_wall) {
        	return(0);
        }
	else {
		return(-(this->arglist_pt->lambda+lambda_t(i,j))*deriv_y_temp(i,j));
	}
}

//divergence du flux de chaleur turbulent (en cartésiennes)
data_t Diff_Eq_Solver::diver_heat_flux_turb(int i, int j) {
	return(deriv_x_heat_flux_x_turb(i,j)+deriv_y_heat_flux_y_turb(i,j));
}

//les dérivées du flux de chaleur turbulent (en cartésiennes)
//attention à ne pas sortir du tableau avec la double dérivée
data_t Diff_Eq_Solver::deriv_x_heat_flux_x_turb(int i, int j) {
	if (is_in(i+2,j) && is_in(i-2,j)) {
		return(
		(heat_flux_x_turb(i+1,j)-heat_flux_x_turb(i-1,j))
		/ (this->arglist_pt->space_step*2)
		);
	}
	else {
		return(0);
	}
}

data_t Diff_Eq_Solver::deriv_y_heat_flux_y_turb(int i, int j) {
	if (is_in(i,j+2) && is_in(i,j-2)) {
		return(
		(heat_flux_y_turb(i,j+1)-heat_flux_y_turb(i,j-1))
		/ (this->arglist_pt->space_step*2)
		);
	}
	else {
	return(0);
	}
}

//la divergence de contrainte_turbulente*vitesse  (cartésiennes)
data_t Diff_Eq_Solver::diver_vtau_turb(int i, int j) {
	return(
	deriv_x_vtaux_turb(i,j)+deriv_y_vtauy_turb(i,j)
	);
}

//définition de contrainte_turbulente*vitesse (cartésiennes)
data_t Diff_Eq_Solver::vtaux_turb(int i, int j) {
	return(
	mesh_grid_1[i][j].speed[0]*tot_stress_xx(i,j)+mesh_grid_1[i][j].speed[1]*tot_stress_xy(i,j)
	);
}

data_t Diff_Eq_Solver::vtauy_turb(int i, int j) {
	return(
	mesh_grid_1[i][j].speed[0]*tot_stress_xy(i,j)+mesh_grid_1[i][j].speed[1]*tot_stress_yy(i,j)
	);
}

//dérivées de la contrainte_turbulente*vitesse (cartésiennes)
//attention à ne pas sortir du tableau
data_t Diff_Eq_Solver::deriv_x_vtaux_turb(int i, int j) {
	if (is_in(i+2,j) && is_in(i-2,j)) {
		return(
		(vtaux_turb(i+1,j)-vtaux_turb(i-1,j))
		/ (this->arglist_pt->space_step*2)
		);
	}
	else {
		return(0);
	}
}

data_t Diff_Eq_Solver::deriv_y_vtauy_turb(int i, int j) {
	if (is_in(i,j+2) && is_in(i,j-2)) {
		return(
		(vtauy_turb(i,j+1)-vtauy_turb(i,j-1))
		/ (this->arglist_pt->space_step*2)
		);
	}
	else {
		return(0);
	}
}

data_t Diff_Eq_Solver::deriv_y_tauxy_turb(int i, int j) {
	if (is_in(i,j+2) && is_in(i,j-2)) {
		return(
		(tot_stress_xy(i,j+1)-tot_stress_xy(i,j-1))
		/ (this->arglist_pt->space_step*2)
		);
	}
	else {
		return(0);
	}
}

data_t Diff_Eq_Solver::deriv_x_tauxx_turb(int i, int j) {
	if (is_in(i+2,j) && is_in(i-2,j)) {
		return(
		(tot_stress_xx(i+1,j)-tot_stress_xx(i-1,j))
		/ (this->arglist_pt->space_step*2)
		);
	}
	else {
		return(0);
	}
}

data_t Diff_Eq_Solver::deriv_y_tauyy_turb(int i, int j) {
	if (is_in(i,j+2) && is_in(i,j-2)) {
		return(
		(tot_stress_yy(i,j+1)-tot_stress_yy(i,j-1))
		/(this->arglist_pt->space_step*2)
		);
	}
	else {
		return(0);
	}
}

data_t Diff_Eq_Solver::deriv_x_tauyx_turb(int i, int j) {
	if (is_in(i+2,j) && is_in(i-2,j)) {
		return(
		(tot_stress_xy(i+1,j)-tot_stress_xy(i-1,j))
		/ (this->arglist_pt->space_step*2)
		);
	}
	else {
		return(0);
	}
}

//définition de la contrainte totale (moléculaire + turbulente)
data_t Diff_Eq_Solver::tot_stress_yy(int i, int j) {
	return(mol_stress_yy(i,j)+turb_stress_yy(i,j));
}

data_t Diff_Eq_Solver::tot_stress_xx(int i, int j) {
	return(mol_stress_xx(i,j)+turb_stress_xx(i,j));
}

data_t Diff_Eq_Solver::tot_stress_xy(int i, int j) {
	return(mol_stress_xy(i,j)+turb_stress_xy(i,j));
}

//définition de la contrainte turbulente (analogie loi de Stokes)
data_t Diff_Eq_Solver::turb_stress_yy(int i, int j) {
	return(2.0*mu_t(i,j)*(2.0/3.0*strain_yy(i,j)-1.0/3.0*strain_xx(i,j))-2.0/3.0*mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].turb_en);
}

data_t Diff_Eq_Solver::turb_stress_xx(int i, int j) {
	return(2.0*mu_t(i,j)*(2.0/3.0*strain_xx(i,j)-1.0/3.0*strain_yy(i,j))-2.0/3.0*mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].turb_en);
}

data_t Diff_Eq_Solver::turb_stress_xy(int i, int j) {
	return(2.0*mu_t(i,j)*strain_xy(i,j));
}

//le "stress" moléculaire : c'est la contrainte moléculaire classique (loi de Stokes)
//attention que la contrainte est symétrique, mol_stress_yx=mol_stress_xy
data_t Diff_Eq_Solver::mol_stress_xy(int i, int j) {
	return(2.0*this->arglist_pt->dyn_visc*strain_xy(i,j));
}

data_t Diff_Eq_Solver::mol_stress_xx(int i, int j) {
	return(2.0*this->arglist_pt->dyn_visc*(2.0/3.0*strain_xx(i,j)-1.0/3.0*strain_yy(i,j)));
}

data_t Diff_Eq_Solver::mol_stress_yy(int i, int j) {
	return(2.0*this->arglist_pt->dyn_visc*(2.0/3.0*strain_yy(i,j)-1.0/3.0*strain_xx(i,j)));
}

//divergence de rho*v*k (modèle k-epsilon) (cartésiennes)
data_t Diff_Eq_Solver::diver_rhovk(int i, int j) {
	return(deriv_x_rhovxk(i,j)+deriv_y_rhovyk(i,j));
}

//dérivées de rho*v_x,y*k (cartésiennes)
data_t Diff_Eq_Solver::deriv_x_rhovxk(int i, int j) {
	return(
	(mesh_grid_1[i+1][j].turb_en*mesh_grid_1[i+1][j].vol_mass*mesh_grid_1[i+1][j].speed[0]-mesh_grid_1[i-1][j].turb_en*mesh_grid_1[i-1][j].vol_mass*mesh_grid_1[i-1][j].speed[0])
	/(2*this->arglist_pt->space_step)
	);
}

data_t Diff_Eq_Solver::deriv_y_rhovyk(int i, int j) {
	return(
	(mesh_grid_1[i][j+1].turb_en*mesh_grid_1[i][j+1].vol_mass*mesh_grid_1[i][j+1].speed[1]-mesh_grid_1[i][j-1].turb_en*mesh_grid_1[i][j-1].vol_mass*mesh_grid_1[i][j-1].speed[1])
	/(2*this->arglist_pt->space_step)
	);
}

//divergence de mu_tot*grad(k) (cartésiennes)
data_t Diff_Eq_Solver::diver_mudk(int i, int j) {
	return(deriv_x_mudxk(i,j)+deriv_y_mudyk(i,j));
}

//dérivées de mu_tot*grad(k) (cartésiennes)
data_t Diff_Eq_Solver::deriv_x_mudxk(int i, int j) {
	return(
	( (mu_t(i+1,j)/sig_k+this->arglist_pt->dyn_visc)*deriv_x_k(i+1,j) - (mu_t(i-1,j)/sig_k+this->arglist_pt->dyn_visc)*deriv_x_k(i-1,j))
	/ (2*this->arglist_pt->space_step)
	);
}

data_t Diff_Eq_Solver::deriv_y_mudyk(int i, int j) {
	return(
	( (mu_t(i,j+1)/sig_k+this->arglist_pt->dyn_visc)*deriv_y_k(i,j+1) - (mu_t(i,j-1)/sig_k+this->arglist_pt->dyn_visc)*deriv_y_k(i,j-1))
	/ (2*this->arglist_pt->space_step)
	);
}

//dérivées de k (cartésiennes)
data_t Diff_Eq_Solver::deriv_x_k(int i, int j) {
	if (mesh_grid_1[i][j].is_wall) {
        	return(0);
        }
	if (is_in(i+1,j) && is_in(i-1,j)) {
		return(
		(mesh_grid_1[i+1][j].turb_en-mesh_grid_1[i-1][j].turb_en)
		/ (2*this->arglist_pt->space_step)
		);
	}
	else {
		return(0.);
	}
}

data_t Diff_Eq_Solver::deriv_y_k(int i, int j) {
	if (mesh_grid_1[i][j].is_wall) {
        	return(0);
        }
	if (is_in(i,j+1) && is_in(i,j-1)) {
		return(
		(mesh_grid_1[i][j+1].turb_en-mesh_grid_1[i][j-1].turb_en)
		/ (2*this->arglist_pt->space_step)
		);
	}
	else {
		return(0.);
	}
}

//divergence de rho*v*epsilon (modèle k-epsilon) (cartésiennes)
data_t Diff_Eq_Solver::diver_rhovepsilon(int i, int j) {
	return(deriv_x_rhovxepsilon(i,j)+deriv_y_rhovyepsilon(i,j));
}

//dérivées de rho*v*epsilon (cartésiennes)
data_t Diff_Eq_Solver::deriv_x_rhovxepsilon(int i, int j) {
	return(
	(mesh_grid_1[i+1][j].turb_dis*mesh_grid_1[i+1][j].vol_mass*mesh_grid_1[i+1][j].speed[0]-mesh_grid_1[i-1][j].turb_dis*mesh_grid_1[i-1][j].vol_mass*mesh_grid_1[i-1][j].speed[0])
	/(2*this->arglist_pt->space_step)
	);
}

data_t Diff_Eq_Solver::deriv_y_rhovyepsilon(int i, int j) {
	return(
	(mesh_grid_1[i][j+1].turb_dis*mesh_grid_1[i][j+1].vol_mass*mesh_grid_1[i][j+1].speed[1]-mesh_grid_1[i][j-1].turb_dis*mesh_grid_1[i][j-1].vol_mass*mesh_grid_1[i][j-1].speed[1])
	/(2*this->arglist_pt->space_step)
	);
}

//divergence de mu_tot*grad(epsilon) (cartésiennes)
data_t Diff_Eq_Solver::diver_mudepsilon(int i, int j) {
	return(deriv_x_mudxepsilon(i,j)+deriv_y_mudyepsilon(i,j));
}

//dérivées de mu_tot*grad(epsilon) (cartésiennes)
data_t Diff_Eq_Solver::deriv_x_mudxepsilon(int i, int j) {
	return(
	((mu_t(i+1,j)/sig_e+this->arglist_pt->dyn_visc)*deriv_x_epsilon(i+1,j) - (mu_t(i-1,j)/sig_e+this->arglist_pt->dyn_visc)*deriv_x_epsilon(i-1,j))
	/ (2*this->arglist_pt->space_step)
	);
}

data_t Diff_Eq_Solver::deriv_y_mudyepsilon(int i, int j) {
	return(
	( (mu_t(i,j+1)/sig_e+this->arglist_pt->dyn_visc)*deriv_y_epsilon(i,j+1) - (mu_t(i,j-1)/sig_e+this->arglist_pt->dyn_visc)*deriv_y_epsilon(i,j-1))
	/ (2*this->arglist_pt->space_step)
	);
}

//dérivées de epsilon (cartésiennes)
data_t Diff_Eq_Solver::deriv_x_epsilon(int i, int j) {
	if (mesh_grid_1[i][j].is_wall) {
        	return(0);
        }
	if (is_in(i+1,j) && is_in(i-1,j)) {
		return(
		(mesh_grid_1[i+1][j].turb_dis-mesh_grid_1[i-1][j].turb_dis)
		/ (2*this->arglist_pt->space_step)
		);
	}
	else {
		return(0.);
	}
}

data_t Diff_Eq_Solver::deriv_y_epsilon(int i, int j) {
	if (mesh_grid_1[i][j].is_wall) {
        	return(0);
        }
	if (is_in(i,j+1) && is_in(i,j-1)) {
		return(
		(mesh_grid_1[i][j+1].turb_dis-mesh_grid_1[i][j-1].turb_dis)
		/ (2*this->arglist_pt->space_step)
		);
	}
	else {
		return(0.);
	}
}

//petite fonction pour savoir si on est dans le tableau ou si on sort.
//si on sort, tous les gradients sont nuls.
bool Diff_Eq_Solver::is_in(int i, int j) {
	if (i>=0 && i<this->arglist_pt->x_size && j>=0 && j<this->arglist_pt->y_size) {
		return(true);
	}
	else {
		return(false);
	}
}

//les update de la masse volumique
//en cartésiennes
void Diff_Eq_Solver::update_vol_mass(int i, int j) {
        mesh_grid_2[i][j].vol_mass = mesh_grid_1[i][j].vol_mass - diver_rhov_c(i,j)*time_step;
}

//en cylindriques
void Diff_Eq_Solver::update_vol_mass_cyl(int i, int j) {
        mesh_grid_2[i][j].vol_mass = mesh_grid_1[i][j].vol_mass-time_step*(1.0/r(i)*deriv_r_rrhovr(i,j)+deriv_y_rhovy(i,j));
}

//les deux updates de la vitesse en cartésiennes
void Diff_Eq_Solver::update_speed_x(int i, int j) {
	mesh_grid_2[i][j].speed[0] = 1.0/mesh_grid_2[i][j].vol_mass*(mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].speed[0]-time_step*(mesh_grid_1[i][j].speed[0]*diver_rhov_c(i,j)+mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].speed[0]*deriv_x_vx(i,j)+mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].speed[1]*deriv_y_vx(i,j)+deriv_x_pres(i,j)-deriv_y_tauxy(i,j)-deriv_x_tauxx(i,j)));
}

void Diff_Eq_Solver::update_speed_y(int i, int j) {
	mesh_grid_2[i][j].speed[1] = 1.0/mesh_grid_2[i][j].vol_mass*(mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].speed[1]-time_step*(mesh_grid_1[i][j].speed[1]*diver_rhov_c(i,j)+mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].speed[0]*deriv_x_vy(i,j)+mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].speed[1]*deriv_y_vy(i,j)+deriv_y_pres(i,j)-deriv_x_tauyx(i,j)-deriv_y_tauyy(i,j)));
}

//les deux updates de la vitesse en cylindriques
void Diff_Eq_Solver::update_speed_r_cyl(int i, int j) {
	mesh_grid_2[i][j].speed[0] = 1.0/mesh_grid_2[i][j].vol_mass*(mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].speed[0]-time_step*(-deriv_x_pres(i,j)-mesh_grid_1[i][j].speed[0]*deriv_x_rhovx(i,j)+mesh_grid_1[i][j].speed[1]*deriv_y_rhovy(i,j)));
}

void Diff_Eq_Solver::update_speed_z_cyl(int i, int j) {
        mesh_grid_2[i][j].speed[1] = 1.0/mesh_grid_2[i][j].vol_mass*(mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].speed[1]-time_step*(deriv_y_pres(i,j) - mesh_grid_1[i][j].speed[0]*deriv_x_rhovx(i,j)+mesh_grid_1[i][j].speed[1]*deriv_y_rhovy(i,j)));
}

//les deux updates de la vitesse en cartésiennes (turbulentes)
void Diff_Eq_Solver::update_speed_x_turb(int i, int j) {
	mesh_grid_2[i][j].speed[0] = 1.0/mesh_grid_2[i][j].vol_mass*(mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].speed[0]-time_step*(mesh_grid_1[i][j].speed[0]*diver_rhov_c(i,j)+mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].speed[0]*deriv_x_vx(i,j)+mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].speed[1]*deriv_y_vx(i,j)+deriv_x_pres(i,j)-deriv_y_tauxy_turb(i,j)-deriv_x_tauxx_turb(i,j)));
}

void Diff_Eq_Solver::update_speed_y_turb(int i, int j) {
	mesh_grid_2[i][j].speed[1] = 1.0/mesh_grid_2[i][j].vol_mass*(mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].speed[1]-time_step*(mesh_grid_1[i][j].speed[1]*diver_rhov_c(i,j)+mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].speed[0]*deriv_x_vy(i,j)+mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].speed[1]*deriv_y_vy(i,j)+deriv_y_pres(i,j)-deriv_x_tauyx_turb(i,j)-deriv_y_tauyy_turb(i,j)));
}

//l'update de la pression pour le gaz de Van der Waals
void Diff_Eq_Solver::update_pres_VDW(int i, int j) {
        mesh_grid_2[i][j].pressure = mesh_grid_2[i][j].vol_mass * R * mesh_grid_2[i][j].temperature / (this->arglist_pt->mol_mass - this->arglist_pt->VDW_b_coef * mesh_grid_2[i][j].vol_mass) - this->arglist_pt->VDW_a_coef * pow(mesh_grid_2[i][j].vol_mass / this->arglist_pt->mol_mass,2);
}

//l'update de la temperature pour le gaz de Van der Waals (cartésiennes)
void Diff_Eq_Solver::update_temp_VDW(int i, int j) {
        mesh_grid_2[i][j].temperature = 2.0/5.0*this->arglist_pt->mol_mass/R*(this->arglist_pt->VDW_a_coef*mesh_grid_2[i][j].vol_mass/pow(this->arglist_pt->mol_mass,2)-1.0/2.0*speed2(i,j,this->mesh_grid_pt1)+1.0/mesh_grid_2[i][j].vol_mass*(mesh_grid_1[i][j].vol_mass*en_tot_VDW(i,j)+time_step*(this->arglist_pt->lambda*(deriv2_x_temp(i,j)+deriv2_y_temp(i,j))-deriv_x_prestotVDWvx(i,j)-deriv_y_prestotVDWvy(i,j))));
}

//l'update de la température pour le gaz de Van der Waals (cylindriques)
void Diff_Eq_Solver::update_temp_VDW_cyl(int i, int j) {
        mesh_grid_2[i][j].temperature = 2.0/5.0*this->arglist_pt->mol_mass/R*(this->arglist_pt->VDW_a_coef*mesh_grid_2[i][j].vol_mass/pow(this->arglist_pt->mol_mass,2)-1.0/2.0*speed2(i,j,this->mesh_grid_pt1)+1.0/mesh_grid_2[i][j].vol_mass*(mesh_grid_1[i][j].vol_mass*en_tot_VDW(i,j)+time_step*(this->arglist_pt->lambda*(-deriv2_x_temp(i,j)+1.0/r(i)*deriv_r_temp(i,j)+deriv2_y_temp(i,j))-1.0/r(i)*deriv_r_prestotVDWrvr(i,j)-deriv_y_prestotVDWvy(i,j))));   
}

//l'update de la température pour le gaz parfait (cartésiennes)
void Diff_Eq_Solver::update_temp_PG(int i, int j) {
	mesh_grid_2[i][j].temperature = 2.0/5.0*this->arglist_pt->mol_mass/R*(-1.0/2.0*speed2(i,j,this->mesh_grid_pt1)+1.0/mesh_grid_2[i][j].vol_mass*(mesh_grid_1[i][j].vol_mass*en_tot_PG(i,j)+time_step*(this->arglist_pt->lambda*(deriv2_x_temp(i,j)+deriv2_y_temp(i,j))-deriv_x_prestotPGvx(i,j)-deriv_y_prestotPGvy(i,j)+diver_vtau(i,j))));
}

//l'update de la température pour le gaz parfait (cylindriques)
/* void Diff_Eq_Solver::update_temp_PG_cyl(int i, int j) {
        mesh_grid_2[i][j].speed[0] = mesh_grid_2[k][l].speed[0];
        mesh_grid_2[i][j].speed[1] = mesh_grid_2[k][l].speed[1];
        mesh_grid_2[i][j].turb_en = mesh_grid_2[k][l].turb_en;
        mesh_grid_2[i][j].turb_dis = mesh_grid_2[k][l].turb_dis;
}*/

//l'update de la température pour le gaz parfait (cylindriques)
void Diff_Eq_Solver::update_temp_PG_cyl(int i, int j) {
	mesh_grid_2[i][j].temperature = 2.0/5.0*this->arglist_pt->mol_mass/R*(-1.0/2.0*speed2(i,j,this->mesh_grid_pt1)+1.0/mesh_grid_2[i][j].vol_mass*(mesh_grid_1[i][j].vol_mass*en_tot_PG(i,j)+time_step*(this->arglist_pt->lambda*(-deriv2_x_temp(i,j)+1.0/r(i)*deriv_r_temp(i,j)+deriv2_y_temp(i,j))-1.0/r(i)*deriv_r_prestotPGrvr(i,j)-deriv_y_prestotPGvy(i,j))));
}
//l'update de la température pour le gaz parfait turbulent (cartésiennes)
void Diff_Eq_Solver::update_temp_PG_turb(int i, int j) { //mise à jour de la température gaz parfait
        mesh_grid_2[i][j].temperature = 2.0/5.0*this->arglist_pt->mol_mass/R*(-1.0/2.0*speed2(i,j,this->mesh_grid_pt1)+1.0/mesh_grid_2[i][j].vol_mass*(mesh_grid_1[i][j].vol_mass*en_tot_PG(i,j)+time_step*(diver_heat_flux_turb(i,j)-deriv_x_prestotPGvx(i,j)-deriv_y_prestotPGvy(i,j)+diver_vtau_turb(i,j))));
}

//l'update de la pression pour le gaz parfait
void Diff_Eq_Solver::update_pres_PG(int i, int j) { //mise à jour de la pression gaz parfait
        mesh_grid_2[i][j].pressure = mesh_grid_2[i][j].vol_mass * R * mesh_grid_2[i][j].temperature / this->arglist_pt->mol_mass;
}

//l'update de l'énergie turbulente (cartésiennes)
void Diff_Eq_Solver::update_k_PG_turb(int i, int j) {
	mesh_grid_2[i][j].turb_en = 1./mesh_grid_2[i][j].vol_mass*(mesh_grid_1[i][j].turb_en*mesh_grid_1[i][j].vol_mass+time_step*(turb_stress_xx(i,j)*strain_xx(i,j)+turb_stress_yy(i,j)*strain_yy(i,j)+2*turb_stress_xy(i,j)*strain_xy(i,j)-mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].turb_dis-diver_rhovk(i,j)+diver_mudk(i,j)));
}

//l'update de la dissipation turbulente (cartésiennes)
void Diff_Eq_Solver::update_epsilon_PG_turb(int i, int j) {

      mesh_grid_2[i][j].turb_dis = 1./mesh_grid_2[i][j].vol_mass*(mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].turb_dis+time_step*(c_e_1*mesh_grid_1[i][j].turb_dis/mesh_grid_1[i][j].turb_en*(turb_stress_xx(i,j)*strain_xx(i,j)+turb_stress_yy(i,j)*strain_yy(i,j)+2*turb_stress_xy(i,j)*strain_xy(i,j)) - c_e_2*(1-0.3*exp(-pow(Ret(i,j),2)))*mesh_grid_1[i][j].vol_mass*pow(mesh_grid_1[i][j].turb_dis,2)/mesh_grid_1[i][j].turb_en - diver_rhovepsilon(i,j) + diver_mudepsilon(i,j)));
}

//petite fonction annexe de copie de ???
void Diff_Eq_Solver::copy_case(int i, int j, int k, int l) {

            mesh_grid_2[i][j].vol_mass = mesh_grid_2[k][l].vol_mass;
            mesh_grid_2[i][j].temperature = mesh_grid_2[k][l].temperature;
            mesh_grid_2[i][j].pressure = mesh_grid_2[k][l].pressure;
            mesh_grid_2[i][j].speed[0] = mesh_grid_2[k][l].speed[0];
            mesh_grid_2[i][j].speed[1] = mesh_grid_2[k][l].speed[1];
            mesh_grid_2[i][j].turb_en = mesh_grid_2[k][l].turb_en;
            mesh_grid_2[i][j].turb_dis = mesh_grid_2[k][l].turb_dis;
}

//petite fonction annexe de calcul de la poussée sur un mesh_grid (intégrale de v_y sur le mur de gauche)
data_t Diff_Eq_Solver::save_thrust() {
  	data_t thrusty=0;
   	for (int k = 0; k<this->arglist_pt->x_size; k++) {
      		thrusty+= -mesh_grid_1[k][0].speed[1]*pow(this->arglist_pt->space_step,2);
   	}
    return(thrusty);
}

//petite fonction annexe qui échange les deux pointeurs vers les deux mesh_grid : permet de passer à l'étape suivante sans perdre l'étape précédente (utile pour le calcul discret)
void Diff_Eq_Solver::exchange_mesh_grid_pts() {
   mesh_grid_t *temp_point = this->mesh_grid_pt2;
   this->mesh_grid_pt2= this->mesh_grid_pt1;
   this->mesh_grid_pt1 = temp_point; //échange des deux pointeurs
}

//Fonctions d'itération

//calcule une itération temporelle des équations différentielles (gaz parfait, cartésiennes)
void Diff_Eq_Solver::calc_iteration_PG_cart() {
  	register int i,j;
  	for (i = 1 ; i < this->arglist_pt->x_size-1; i++) {
    		for (j = 1 ; j < this->arglist_pt->y_size-1; j++) {
			if (not(mesh_grid_1[i][j].is_wall)) {
        		this->update_vol_mass(i,j);
        		this->update_speed_x(i,j);
        		this->update_speed_y(i,j);
        		this->update_temp_PG(i,j);
        		this->update_pres_PG(i,j);
			}
      		}
   	}
   	for (i = 1 ; i < this->arglist_pt->x_size-1; i++) {
        	copy_case(i,0,i,1);
   	}
   	for (j = 1 ; j < this->arglist_pt->y_size-1; j++) {
           copy_case(0,j,1,j);
           copy_case(this->arglist_pt->x_size-1,j,this->arglist_pt->x_size-2,j);
   	}
   	this->exchange_mesh_grid_pts();
}
















//calcule une itération temporelle des équations différentielles (gaz parfait, cartésiennes, turbulentes)

void Diff_Eq_Solver::partial_calc_iteration_PG_cart_turb(int i_min, int i_max, int thread_id)
{
    register int i,j;
    for (i = i_min ; i < i_max; i++) {
        for (j = 1 ; j < this->arglist_pt->y_size-1; j++) {
            if (not(mesh_grid_1[i][j].is_wall)) {
                this->update_vol_mass(i,j);
                this->update_speed_x_turb(i,j);
                this->update_speed_y_turb(i,j);
                this->update_temp_PG_turb(i,j);
                this->update_pres_PG(i,j);
                this->update_k_PG_turb(i,j);
                this->update_epsilon_PG_turb(i,j);
            }
        }
    }
}





void Diff_Eq_Solver::calc_iteration_PG_cart_turb() {
    // on découpe la matrice en bandes horizontales -> chaque calcul part sur un coeur
    register int i, j, t, i_min, i_max, slice;
    i_min = 1;
    slice = (this->arglist_pt->x_size - 1)/ this->arglist_pt->nb_of_threads; // taille d'une tranche horizontale
    i_max = i_min + slice;

    for (t = 0 ; t < this->arglist_pt->nb_of_threads - 1; t++){ 
        this->threads[t] = thread (&Diff_Eq_Solver::partial_calc_iteration_PG_cart_turb, this, i_min, i_max, t); // on envoie sur les coeurs
        i_min = i_max;
        i_max += slice;
    }

    this->threads[this->arglist_pt->nb_of_threads - 1] = 
        thread (&Diff_Eq_Solver::partial_calc_iteration_PG_cart_turb, this, i_min, this->arglist_pt->x_size-1, this->arglist_pt->nb_of_threads-1); // le dernier se tape le rab de la division euclidienne

    for (auto &thrd : this->threads) thrd.join(); 

    
    for (i = 1 ; i < this->arglist_pt->x_size-1; i++) {
           	copy_case(i,0,i,1);
   	}
   	for (j = 1 ; j < this->arglist_pt->y_size-1; j++) {
           	copy_case(0,j,1,j);
           	copy_case(this->arglist_pt->x_size-1,j,this->arglist_pt->x_size-2,j);
   	}
}



















//calcule une itération temporelle des équations différentielles (gaz parfait, cylindriques)
void Diff_Eq_Solver::calc_iteration_PG_cyl() {
        register int i,j;
        for (i = 1 ; i < this->arglist_pt->x_size-1; i++) {
                for (j = 1 ; j < this->arglist_pt->y_size-1; j++) {
      			if (not(mesh_grid_1[i][j].is_wall)) {
                        this->update_vol_mass_cyl(i,j);
                        this->update_speed_r_cyl(i,j);
                        this->update_speed_z_cyl(i,j);
                        this->update_temp_PG_cyl(i,j);
                        this->update_pres_PG(i,j);
      			}
                }
        }
        for (i = 1 ; i < this->arglist_pt->x_size-1; i++) {
                copy_case(i,0,i,1);
        }
        for (j = 1 ; j < this->arglist_pt->y_size-1; j++) {
                copy_case(0,j,1,j);
                copy_case(this->arglist_pt->x_size-1,j,this->arglist_pt->x_size-2,j);
        }
        this->exchange_mesh_grid_pts();
}

//calcule une itération temporelle des équations différentielles (gaz de Van der Waals, cartésiennes)
void Diff_Eq_Solver::calc_iteration_VDW_cart() {
        register int i,j;
        for (i = 1 ; i < this->arglist_pt->x_size-1; i++) {
                for (j = 1 ; j < this->arglist_pt->y_size-1; j++) {
      			if (not(mesh_grid_1[i][j].is_wall)) {
                        this->update_vol_mass(i,j);
                        this->update_speed_x(i,j);
                        this->update_speed_y(i,j);
                        this->update_temp_VDW(i,j);
                        this->update_pres_VDW(i,j);
      			}
                }
        }
        for (i = 1 ; i < this->arglist_pt->x_size-1; i++) {
                copy_case(i,0,i,1);
        }
        for (j = 1 ; j < this->arglist_pt->y_size-1; j++) {
                copy_case(0,j,1,j);
                copy_case(this->arglist_pt->x_size-1,j,this->arglist_pt->x_size-2,j);
        }
        this->exchange_mesh_grid_pts();
}

//calcule une itération temporelle des équations différentielles (gaz de Van der Waals, cylindriques)
void Diff_Eq_Solver::calc_iteration_VDW_cyl() {
        register int i,j;
        for (i = 1 ; i < this->arglist_pt->x_size-1; i++) {
                for (j = 1 ; j < this->arglist_pt->y_size-1; j++) {
      			if (not(mesh_grid_1[i][j].is_wall)) {
                        this->update_vol_mass_cyl(i,j);
                        this->update_speed_r_cyl(i,j);
                        this->update_speed_z_cyl(i,j);
                        this->update_temp_VDW_cyl(i,j);
                        this->update_pres_VDW(i,j);
      			}
                }
        }
        for (i = 1 ; i < this->arglist_pt->x_size-1; i++) {
                copy_case(i,0,i,1);
        }
        for (j = 1 ; j < this->arglist_pt->y_size-1; j++) {
                copy_case(0,j,1,j);
                copy_case(this->arglist_pt->x_size-1,j,this->arglist_pt->x_size-2,j);
        }
        this->exchange_mesh_grid_pts();
}

//lance l'itération de toutes les étapes du diff_eq_solver (gaz parfait, cartésiennes)
void Diff_Eq_Solver::solve_PG_cart() {
        register int i;
        for (i=0; i<this->arglist_pt->iter_number_solver; i++) {
        	this->calc_iteration_PG_cart();
   		this->thrust[i] = save_thrust();
        }
        this->DM->temporal_stuff_plotter(this);
}

//lance l'itération de toutes les étapes du diff_eq_solver (gaz parfait, cartésiennes, turbulent, gradient initial)
void Diff_Eq_Solver::solve_PG_cart_turb_init_grad() {
        register double t_step;
        this->UI->start_DES_chrono();

        for (this->ite_count=0; this->ite_count < this->arglist_pt->iter_number_solver; this->ite_count++) {
                this->calc_iteration_PG_cart_turb();
                this->find_all_max();
                this->calc_all_max(this->ite_count);
                t_step = sqrt( this->arglist_pt->space_step /
                            ( (pow(this->variables_max.speed0[this->ite_count],2) 
                                + pow(this->variables_max.speed1[this->ite_count], 2) ) 
                            * this->arglist_pt->CFL_cond) );
                if ( t_step > this->arglist_pt->init_time_step ) t_step = this->arglist_pt->init_time_step;
                this->time_steps[this->ite_count] = t_step;  
                this->thrust[this->ite_count] = save_thrust();
                this->exchange_mesh_grid_pts();
                this->UI->refresh_DES_progress(this->ite_count, this->arglist_pt->iter_number_solver);
        }
        this->ite_count--;        
        this->UI->space();
        this->DM->temporal_stuff_plotter(this);
}

//lance l'itération de toutes les étapes du diff_eq_solver (gaz parfait, cartésiennes, turbulent, conditions au fond de la chambre en évolution)
void Diff_Eq_Solver::solve_PG_cart_turb_evol_chamber() {
        register double t_step;
        this->UI->start_DES_chrono();

        for (this->ite_count=0; this->ite_count < this->arglist_pt->iter_number_evol_chamber; this->ite_count++) {
                this->calc_iteration_PG_cart_turb();
                this->NP->update_chamber_cond();
                this->find_all_max();
                this->calc_all_max(this->ite_count);
                t_step = sqrt( this->arglist_pt->space_step /
                            ( (pow(this->variables_max.speed0[this->ite_count],2) 
                                + pow(this->variables_max.speed1[this->ite_count], 2) ) 
                            * this->arglist_pt->CFL_cond) );
                if ( t_step > this->arglist_pt->init_time_step ) t_step = this->arglist_pt->init_time_step;
                this->time_steps[this->ite_count] = t_step;  
                this->thrust[this->ite_count] = save_thrust();
                this->exchange_mesh_grid_pts();
                this->UI->refresh_DES_progress(this->ite_count, this->arglist_pt->iter_number_solver);
        }

        for (this->ite_count = this->ite_count; this->ite_count < this->arglist_pt->iter_number_solver; this->ite_count++) {
                this->calc_iteration_PG_cart_turb();
                this->find_all_max();
                this->calc_all_max(this->ite_count);
                t_step = sqrt( this->arglist_pt->space_step /
                            ( (pow(this->variables_max.speed0[this->ite_count],2) 
                                + pow(this->variables_max.speed1[this->ite_count], 2) ) 
                            * this->arglist_pt->CFL_cond) );
                if ( t_step > this->arglist_pt->init_time_step ) t_step = this->arglist_pt->init_time_step;
                this->time_steps[this->ite_count] = t_step;  
                this->thrust[this->ite_count] = save_thrust();
                this->exchange_mesh_grid_pts();
                this->UI->refresh_DES_progress(this->ite_count, this->arglist_pt->iter_number_solver);
        }

        this->ite_count--;        
        this->UI->space();
        this->DM->temporal_stuff_plotter(this);
}
//lance l'itération de toutes les étapes du diff_eq_solver (gaz de Van der Waals, cartésiennes)
void Diff_Eq_Solver::solve_VDW_cart() {
        register int i;
        for (i=0; i<this->arglist_pt->iter_number_solver; i++) {
        	this->calc_iteration_VDW_cart();
        	this->thrust[i] = save_thrust();
        }
        this->DM->temporal_stuff_plotter(this);
}

//lance l'itération de toutes les étapes du diff_eq_solver (gaz parfait, cylindriques)
void Diff_Eq_Solver::solve_PG_cyl() {
        register int i;
        for (i=0; i<this->arglist_pt->iter_number_solver; i++) {
        	this->calc_iteration_PG_cyl();
        	this->thrust[i] = save_thrust();
        }
        this->DM->temporal_stuff_plotter(this);
}

//lance l'itération de toutes les étapes du diff_eq_solver (gaz de Van der Waals, cylindriques)
void Diff_Eq_Solver::solve_VDW_cyl() {
        register int i;
        for (i=0; i<this->arglist_pt->iter_number_solver; i++) {
        	this->calc_iteration_VDW_cyl();
        	this->thrust[i] = save_thrust();
        }
        this->DM->temporal_stuff_plotter(this);
}

//Lance le bon solveur selon l'argument diff_eq_solver_algo de l'argfile
void Diff_Eq_Solver::solve()
{
        switch ( this->arglist_pt->diff_eq_solver_algo )
        {
                case PG_cart:
                        this->solve_PG_cart();
                        break;

                case VDW_cart:
                        this->solve_VDW_cart();
                        break;

                case PG_cyl:
                        this->solve_PG_cart();
                        break;

                case VDW_cyl:
                        this->solve_VDW_cyl();
                        break;
                
                case PG_cart_turb:
                        switch (this->arglist_pt->init_cond_type)
                        {
                            case INIT_GRAD:
                                this->solve_PG_cart_turb_init_grad();
                                break;
                            case EVOL_CHAMBER:
                                this->solve_PG_cart_turb_evol_chamber();
                                break;
                        }
                        break;

                default:
                        throw "DES : Bad solving algorithm";
        }
}
