/* This is differential equation solver implementation linked to diff_eq_solver.h
 *
 *
 *
 *
*/


// ======== includes ======== 
#include "../include/diff_eq_solver.h"
#include "../include/arg_interpreter.h"
#include "../include/nozzle_profiler.h"
#include <cmath>
#define R 8.314
#include <fstream>
#include <iostream>

// ======== implementation ========

Diff_Eq_Solver::Diff_Eq_Solver(Usr_Interface *UI, arglist_struct *arglist_pt, mesh_grid_t *mesh_grid_pt1, mesh_grid_t *mesh_grid_pt2) {
        this->UI = UI;
        this->arglist_pt = arglist_pt;
        this->mesh_grid_pt1 = mesh_grid_pt1;
        this->mesh_grid_pt2 = mesh_grid_pt2;
}

data_t Diff_Eq_Solver::speed2(int i, int j, int k) {
        mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
        mesh_grid_t &mesh_grid_2 = (*(this->mesh_grid_pt2));
        
        if (k==1) {
                return(pow(mesh_grid_1[i][j].speed[0],2) + pow(mesh_grid_1[i][j].speed[1],2));
        }
        else {
                return(pow(mesh_grid_2[i][j].speed[0],2) + pow(mesh_grid_2[i][j].speed[1],2));
        }
}

data_t Diff_Eq_Solver::pres_tot_PG(int i, int j) {
        mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
        
        return(
        mesh_grid_1[i][j].vol_mass*en_tot_PG(i,j)+mesh_grid_1[i][j].pressure
        );
}

data_t Diff_Eq_Solver::pres_tot_VDW(int i, int j) {
        mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
        
        return(
        mesh_grid_1[i][j].vol_mass*en_tot_VDW(i,j)+mesh_grid_1[i][j].pressure
        );
}

data_t Diff_Eq_Solver::en_tot_PG(int i, int j) {
        mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
        
        return(
        5.0/2.0*R/this->arglist_pt->mol_mass*mesh_grid_1[i][j].temperature + 1.0/2.0*speed2(i,j,1)
        );
}

data_t Diff_Eq_Solver::en_tot_VDW(int i, int j) {
        mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
        
        return(
        5.0/2.0*R*mesh_grid_1[i][j].temperature/this->arglist_pt->mol_mass + 1.0/2.0*speed2(i,j,1) - this->arglist_pt->VDW_a_coef*mesh_grid_1[i][j].vol_mass/pow(this->arglist_pt->mol_mass,2)
        );
}

data_t Diff_Eq_Solver::diver_rhov_c(int i, int j) {
        mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
        
        return(
             deriv_x_rhovx(i,j)+deriv_y_rhovy(i,j)
        );
}

data_t Diff_Eq_Solver::deriv_x_temp(int i, int j) {
        mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
        
        return(
        (1-mesh_grid_1[i+1][j].is_wall)*(1-mesh_grid_1[i-1][j].is_wall)*(mesh_grid_1[i+1][j].temperature-mesh_grid_1[i-1][j].temperature)
        / (2*this->arglist_pt->space_step)
        );
}

data_t Diff_Eq_Solver::deriv_y_temp(int i, int j) {
        mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
        
        return(
        (1-mesh_grid_1[i][j+1].is_wall)*(1-mesh_grid_1[i][j-1].is_wall)*(mesh_grid_1[i][j+1].temperature-mesh_grid_1[i][j-1].temperature)
        / (2*this->arglist_pt->space_step)
        );
}

data_t Diff_Eq_Solver::deriv2_x_temp(int i, int j) {
        mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
        
        return(
        ((1-mesh_grid_1[i+1][j].is_wall)*(mesh_grid_1[i+1][j].temperature-mesh_grid_1[i][j].temperature)-(1-mesh_grid_1[i-1][j].is_wall)*(mesh_grid_1[i][j].temperature-mesh_grid_1[i-1][j].temperature))
        / (pow(this->arglist_pt->space_step,2))
        );
}

data_t Diff_Eq_Solver::deriv2_y_temp(int i, int j) {
        mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
        
        return(
        ((1-mesh_grid_1[i][j+1].is_wall)*(mesh_grid_1[i][j+1].temperature-mesh_grid_1[i][j].temperature)-(1-mesh_grid_1[i][j-1].is_wall)*(mesh_grid_1[i][j].temperature-mesh_grid_1[i][j-1].temperature))
        / (pow(this->arglist_pt->space_step,2))
        );
}
data_t Diff_Eq_Solver::deriv_x_pres(int i, int j) {
        mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
        
        return(
        (1-mesh_grid_1[i+1][j].is_wall)*(1-mesh_grid_1[i-1][j].is_wall)*(mesh_grid_1[i+1][j].pressure-mesh_grid_1[i-1][j].pressure)
        / (this->arglist_pt->space_step*2)
        );
}

data_t Diff_Eq_Solver::deriv_y_pres(int i, int j) {
        mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
        
        return(
        (1-mesh_grid_1[i][j+1].is_wall)*(1-mesh_grid_1[i][j-1].is_wall)*(mesh_grid_1[i][j+1].pressure-mesh_grid_1[i][j-1].pressure)
        / (this->arglist_pt->space_step*2)
        );
}

data_t Diff_Eq_Solver::deriv_x_rhovx(int i, int j) {
        mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
        
        return(
        (mesh_grid_1[i+1][j].vol_mass*mesh_grid_1[i+1][j].speed[0]-mesh_grid_1[i-1][j].vol_mass*mesh_grid_1[i-1][j].speed[0])
        / (2*this->arglist_pt->space_step)
        );
}

data_t Diff_Eq_Solver::deriv_y_rhovx(int i, int j) {
        mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
        
        return(
        (mesh_grid_1[i][j+1].vol_mass*mesh_grid_1[i][j+1].speed[0]-mesh_grid_1[i][j-1].vol_mass*mesh_grid_1[i][j-1].speed[0])
        / (this->arglist_pt->space_step*2)
        );
}

data_t Diff_Eq_Solver::deriv_x_rhovy(int i, int j) {
        mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
        
        return(
        (mesh_grid_1[i+1][j].vol_mass*mesh_grid_1[i+1][j].speed[1]-mesh_grid_1[i-1][j].vol_mass*mesh_grid_1[i-1][j].speed[1])
        / (this->arglist_pt->space_step*2)
        );
}

data_t Diff_Eq_Solver::deriv_y_rhovy(int i, int j) {
        mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
        
        return(
        (mesh_grid_1[i][j+1].vol_mass*mesh_grid_1[i][j+1].speed[1]-mesh_grid_1[i][j-1].vol_mass*mesh_grid_1[i][j-1].speed[1])
        / (this->arglist_pt->space_step*2)
        );
}

data_t Diff_Eq_Solver::deriv_x_prestotPGvx(int i, int j) {
        mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
        
        return(
        (mesh_grid_1[i+1][j].speed[0]*pres_tot_PG(i+1,j)-mesh_grid_1[i-1][j].speed[0]*pres_tot_PG(i-1,j))
        / (this->arglist_pt->space_step*2)
        );
}

data_t Diff_Eq_Solver::deriv_y_prestotPGvy(int i, int j) {
        mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
        
        return(
        (mesh_grid_1[i][j+1].speed[1]*pres_tot_PG(i,j+1)-mesh_grid_1[i][j-1].speed[1]*pres_tot_PG(i,j-1))
        / (this->arglist_pt->space_step*2)
        );
}

data_t Diff_Eq_Solver::deriv_r_prestotPGrvr(int i, int j) {
        mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
        
        return(
        -(r(i+1)*mesh_grid_1[i+1][j].speed[0]*pres_tot_PG(i+1,j)-r(i-1)*mesh_grid_1[i-1][j].speed[0]*pres_tot_PG(i-1,j))
        / (this->arglist_pt->space_step*2)
        );
}

data_t Diff_Eq_Solver::deriv_r_temp(int i, int j) {
        mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
        
        return(
        (1-mesh_grid_1[i+1][j].is_wall)*(1-mesh_grid_1[i-1][j].is_wall)*(mesh_grid_1[i-1][j].temperature-mesh_grid_1[i+1][j].temperature)
        / (this->arglist_pt->space_step*2)
        );
}

data_t Diff_Eq_Solver::deriv_x_prestotVDWvx(int i, int j) {
        mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
        
        return(
        (mesh_grid_1[i+1][j].speed[0]*pres_tot_VDW(i+1,j)-mesh_grid_1[i-1][j].speed[0]*pres_tot_VDW(i-1,j))
        / (this->arglist_pt->space_step*2)
        );
}

data_t Diff_Eq_Solver::deriv_y_prestotVDWvy(int i, int j) {
        mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));

        return(
        (mesh_grid_1[i][j+1].speed[1]*pres_tot_VDW(i,j+1)-mesh_grid_1[i][j-1].speed[1]*pres_tot_VDW(i,j-1))
        / (this->arglist_pt->space_step*2)
        );
}

data_t Diff_Eq_Solver::deriv_r_prestotVDWrvr(int i, int j) {
        mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
        
        return(
        -(r(i+1)*mesh_grid_1[i+1][j].speed[0]*pres_tot_VDW(i+1,j)-r(i-1)*mesh_grid_1[i-1][j].speed[0]*pres_tot_VDW(i-1,j))
        / (this->arglist_pt->space_step*2)
        );
}

data_t Diff_Eq_Solver::r(int i) {
        return(this->arglist_pt->space_step*(this->arglist_pt->x_size-i));
}

data_t Diff_Eq_Solver::deriv_r_rrhovr(int i, int j) {
        mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
        
        return(
        (r(i-1)*mesh_grid_1[i-1][j].vol_mass*mesh_grid_1[i-1][j].speed[0]-r(i+1)*mesh_grid_1[i+1][j].vol_mass*mesh_grid_1[i+1][j].speed[0])
        / (this->arglist_pt->space_step*2)
        );
}

data_t Diff_Eq_Solver::deriv_x_vx(int i, int j) {
	mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
	
	return(
	(mesh_grid_1[i+1][j].speed[0]-mesh_grid_1[i-1][j].speed[0])
	/ (this->arglist_pt->space_step*2)
	);
}

data_t Diff_Eq_Solver::deriv_y_vx(int i, int j) {
	mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));

	return(
	(mesh_grid_1[i][j+1].speed[0]-mesh_grid_1[i][j-1].speed[0])
	/ (this->arglist_pt->space_step*2)
	);
}

data_t Diff_Eq_Solver::deriv_x_vy(int i, int j) {
	mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));

	return(
	(mesh_grid_1[i+1][j].speed[1]-mesh_grid_1[i-1][j].speed[1])
	/ (this->arglist_pt->space_step*2)
	);
}

data_t Diff_Eq_Solver::deriv_y_vy(int i, int j) {
	mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));

	return(
	(mesh_grid_1[i][j+1].speed[1]-mesh_grid_1[i][j-1].speed[1])
	/ (this->arglist_pt->space_step*2)
	);
}

data_t Diff_Eq_Solver::mol_stress_xy(int i, int j) {
	return(2.0*this->arglist_pt->dyn_visc*strain_xy(i,j));
}

data_t Diff_Eq_Solver::mol_stress_xx(int i, int j) {
	return(2.0*this->arglist_pt->dyn_visc*(2.0/3.0*strain_xx(i,j)-1.0/3.0*strain_yy(i,j)));
}

data_t Diff_Eq_Solver::mol_stress_yy(int i, int j) {
	return(2.0*this->arglist_pt->dyn_visc*(2.0/3.0*strain_yy(i,j)-1.0/3.0*strain_xx(i,j)));
}

data_t Diff_Eq_Solver::strain_xy(int i, int j) {
	return(1.0/2.0*(deriv_y_vx(i,j)+deriv_x_vy(i,j)));
}

data_t Diff_Eq_Solver::strain_xx(int i, int j) {
	return(deriv_x_vx(i,j));
}

data_t Diff_Eq_Solver::strain_yy(int i, int j) {
	return(deriv_y_vy(i,j));
}

data_t Diff_Eq_Solver::deriv_y_tauxy(int i, int j) {
	mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
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
	mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
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
	mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
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
	mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
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

data_t Diff_Eq_Solver::diver_vtau(int i, int j) {
	mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
	return(
	deriv_x_vtaux(i,j)+deriv_y_vtauy(i,j)
	);
}

data_t Diff_Eq_Solver::vtaux(int i, int j) {
	mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
	return(
	mesh_grid_1[i][j].speed[0]*mol_stress_xx(i,j)+mesh_grid_1[i][j].speed[1]*mol_stress_xy(i,j)
	);
}

data_t Diff_Eq_Solver::vtauy(int i, int j) {
	mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
	return(
	mesh_grid_1[i][j].speed[0]*mol_stress_xy(i,j)+mesh_grid_1[i][j].speed[1]*mol_stress_yy(i,j)
	);
}

data_t Diff_Eq_Solver::deriv_x_vtaux(int i, int j) {
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

bool Diff_Eq_Solver::is_in(int i, int j) {
	if (i>=0 && i<this->arglist_pt->x_size && j>=0 && j<this->arglist_pt->y_size) {
		return(true);
	}
	else {
		return(false);
	}
}

void Diff_Eq_Solver::update_vol_mass(int i, int j) { //mise à jour de la masse volumique
        mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
        mesh_grid_t &mesh_grid_2 = (*(this->mesh_grid_pt2));
        
        if (not (mesh_grid_1[i][j].is_wall)) {
                mesh_grid_2[i][j].vol_mass = mesh_grid_1[i][j].vol_mass - diver_rhov_c(i,j)*this->arglist_pt->time_step;
        }
}

void Diff_Eq_Solver::update_vol_mass_cyl(int i, int j) {
        mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
        mesh_grid_t &mesh_grid_2 = (*(this->mesh_grid_pt2));
        
        if (not (mesh_grid_1[i][j].is_wall)){
                mesh_grid_2[i][j].vol_mass = mesh_grid_1[i][j].vol_mass-this->arglist_pt->time_step*(1.0/r(i)*deriv_r_rrhovr(i,j)+deriv_y_rhovy(i,j));
        }
}

void Diff_Eq_Solver::update_speed_x(int i, int j) {
	mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
	mesh_grid_t &mesh_grid_2 = (*(this->mesh_grid_pt2));
	
	if (not(mesh_grid_1[i][j].is_wall)) {
	mesh_grid_2[i][j].speed[0] = 1.0/mesh_grid_2[i][j].vol_mass*(mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].speed[0]-this->arglist_pt->time_step*(mesh_grid_1[i][j].speed[0]*diver_rhov_c(i,j)+mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].speed[0]*deriv_x_vx(i,j)+mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].speed[1]*deriv_y_vx(i,j)+deriv_x_pres(i,j)-deriv_y_tauxy(i,j)-deriv_x_tauxx(i,j) ) );
	}
}

void Diff_Eq_Solver::update_speed_y(int i, int j) {
	mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
	mesh_grid_t &mesh_grid_2 = (*(this->mesh_grid_pt2));

	if (not(mesh_grid_1[i][j].is_wall)) {
	mesh_grid_2[i][j].speed[1] = 1.0/mesh_grid_2[i][j].vol_mass*(mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].speed[1]-this->arglist_pt->time_step*(mesh_grid_1[i][j].speed[1]*diver_rhov_c(i,j)+mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].speed[0]*deriv_x_vy(i,j)+mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].speed[1]*deriv_y_vy(i,j)+deriv_y_pres(i,j)-deriv_x_tauyx(i,j)-deriv_y_tauyy(i,j)  )  );
	}
}

void Diff_Eq_Solver::update_speed_r_cyl(int i, int j) {
        mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
        mesh_grid_t &mesh_grid_2 = (*(this->mesh_grid_pt2));
        
        if (not (mesh_grid_1[i][j].is_wall) ){
                mesh_grid_2[i][j].speed[0] = 1.0/mesh_grid_2[i][j].vol_mass*(mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].speed[0]-this->arglist_pt->time_step*(-deriv_x_pres(i,j)-mesh_grid_1[i][j].speed[0]*deriv_x_rhovx(i,j)+mesh_grid_1[i][j].speed[1]*deriv_y_rhovy(i,j)));
        }
}

void Diff_Eq_Solver::update_speed_z_cyl(int i, int j) {
        mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
        mesh_grid_t &mesh_grid_2 = (*(this->mesh_grid_pt2));
        
        if (not (mesh_grid_1[i][j].is_wall)) {
                mesh_grid_2[i][j].speed[1] = 1.0/mesh_grid_2[i][j].vol_mass*(mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].speed[1]-this->arglist_pt->time_step*(deriv_y_pres(i,j) - mesh_grid_1[i][j].speed[0]*deriv_x_rhovx(i,j)+mesh_grid_1[i][j].speed[1]*deriv_y_rhovy(i,j) )  );
        }
}

void Diff_Eq_Solver::update_pres_VDW(int i, int j) { //mise à jour de la pression Van der Waals
        mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
        mesh_grid_t &mesh_grid_2 = (*(this->mesh_grid_pt2));

        if (not (mesh_grid_1[i][j].is_wall)) {
        mesh_grid_2[i][j].pressure = mesh_grid_2[i][j].vol_mass * R * mesh_grid_2[i][j].temperature / (this->arglist_pt->mol_mass - this->arglist_pt->VDW_b_coef * mesh_grid_2[i][j].vol_mass) - this->arglist_pt->VDW_a_coef * pow(mesh_grid_2[i][j].vol_mass / this->arglist_pt->mol_mass,2);
        }
}

void Diff_Eq_Solver::update_temp_VDW(int i, int j) { //mise à jour de la température Van der Waals
        mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
        mesh_grid_t &mesh_grid_2 = (*(this->mesh_grid_pt2));

        if (not (mesh_grid_1[i][j].is_wall)) {
                mesh_grid_2[i][j].temperature = 2.0/5.0*this->arglist_pt->mol_mass/R*(this->arglist_pt->VDW_a_coef*mesh_grid_2[i][j].vol_mass/pow(this->arglist_pt->mol_mass,2)-1.0/2.0*speed2(i,j,2)+1.0/mesh_grid_2[i][j].vol_mass*(mesh_grid_1[i][j].vol_mass*en_tot_VDW(i,j)+this->arglist_pt->time_step*(this->arglist_pt->lambda*(deriv2_x_temp(i,j)+deriv2_y_temp(i,j))-deriv_x_prestotVDWvx(i,j)-deriv_y_prestotVDWvy(i,j))));
        }
}

void Diff_Eq_Solver::update_temp_VDW_cyl(int i, int j) {
        mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
        mesh_grid_t &mesh_grid_2 = (*(this->mesh_grid_pt2));
        
        if (not (mesh_grid_1[i][j].is_wall) ){
                mesh_grid_2[i][j].temperature = 2.0/5.0*this->arglist_pt->mol_mass/R*(this->arglist_pt->VDW_a_coef*mesh_grid_2[i][j].vol_mass/pow(this->arglist_pt->mol_mass,2)-1.0/2.0*speed2(i,j,2)+1.0/mesh_grid_2[i][j].vol_mass*(mesh_grid_1[i][j].vol_mass*en_tot_VDW(i,j)+this->arglist_pt->time_step*(this->arglist_pt->lambda*(-deriv2_x_temp(i,j)+1.0/r(i)*deriv_r_temp(i,j)+deriv2_y_temp(i,j))-1.0/r(i)*deriv_r_prestotVDWrvr(i,j)-deriv_y_prestotVDWvy(i,j))));   
        }
}

void Diff_Eq_Solver::update_temp_PG(int i, int j) { //mise à jour de la température gaz parfait
        mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
        mesh_grid_t &mesh_grid_2 = (*(this->mesh_grid_pt2));

        if (not (mesh_grid_1[i][j].is_wall) ){
                mesh_grid_2[i][j].temperature = 2.0/5.0*this->arglist_pt->mol_mass/R*(-1.0/2.0*speed2(i,j,2)+1.0/mesh_grid_2[i][j].vol_mass*(mesh_grid_1[i][j].vol_mass*en_tot_PG(i,j)+this->arglist_pt->time_step*(this->arglist_pt->lambda*(deriv2_x_temp(i,j)+deriv2_y_temp(i,j))-deriv_x_prestotPGvx(i,j)-deriv_y_prestotPGvy(i,j)+diver_vtau(i,j))));
        }
}

void Diff_Eq_Solver::update_temp_PG_cyl(int i, int j) {
        mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
        mesh_grid_t &mesh_grid_2 = (*(this->mesh_grid_pt2));
        
        if (not (mesh_grid_1[i][j].is_wall)) {
                mesh_grid_2[i][j].temperature = 2.0/5.0*this->arglist_pt->mol_mass/R*(-1.0/2.0*speed2(i,j,2)+1.0/mesh_grid_2[i][j].vol_mass*(mesh_grid_1[i][j].vol_mass*en_tot_PG(i,j)+this->arglist_pt->time_step*(this->arglist_pt->lambda*(-deriv2_x_temp(i,j)+1.0/r(i)*deriv_r_temp(i,j)+deriv2_y_temp(i,j))-1.0/r(i)*deriv_r_prestotPGrvr(i,j)-deriv_y_prestotPGvy(i,j))));
        }
}

void Diff_Eq_Solver::update_pres_PG(int i, int j) { //mise à jour de la pression gaz parfait
        mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
        mesh_grid_t &mesh_grid_2 = (*(this->mesh_grid_pt2));
        
        if (not (mesh_grid_1[i][j].is_wall)) {
                mesh_grid_2[i][j].pressure = mesh_grid_2[i][j].vol_mass * R * mesh_grid_2[i][j].temperature / this->arglist_pt->mol_mass;
        }
}

void Diff_Eq_Solver::thrust_saver() {
	mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
	
	data_t thrust = 0.;
	for (int i =0;i<this->arglist_pt->x_size;i++) {
		thrust+=-mesh_grid_1[i][0].speed[0]*pow(this->arglist_pt->space_step,2);
	}
	ofstream file;
	file.open("./log.txt",ios::out);
	file << thrust << endl;
	file.close();
}

void Diff_Eq_Solver::copy_case(int i, int j, int k, int l) {
        mesh_grid_t &mesh_grid_2 = (*(this->mesh_grid_pt2));
        
        mesh_grid_2[i][j].vol_mass = mesh_grid_2[k][l].vol_mass;
        mesh_grid_2[i][j].temperature = mesh_grid_2[k][l].temperature;
        mesh_grid_2[i][j].pressure = mesh_grid_2[k][l].pressure;
        mesh_grid_2[i][j].speed[0] = mesh_grid_2[k][l].speed[0];
        mesh_grid_2[i][j].speed[1] = mesh_grid_2[k][l].speed[1];
}




// Une iteration totale

void Diff_Eq_Solver::calc_iteration_PG_cart() {
        mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
        mesh_grid_t &mesh_grid_2 = (*(this->mesh_grid_pt2));
        
  register int i,j;
  for (i = 1 ; i < this->arglist_pt->x_size-1; i++) {
    for (j = 1 ; j < this->arglist_pt->y_size-1; j++) {
      
        this->update_vol_mass(i,j);
        
       this->update_speed_x(i,j);
        
        this->update_speed_y(i,j);
        
        this->update_temp_PG(i,j);
        
        this->update_pres_PG(i,j);
        
      }
   }
   for (i = 1 ; i < this->arglist_pt->x_size-1; i++) {
           copy_case(i,0,i,1);
   }
   for (j = 1 ; j < this->arglist_pt->y_size-1; j++) {
           copy_case(0,j,1,j);
           copy_case(this->arglist_pt->x_size-1,j,this->arglist_pt->x_size-2,j);
   }
   
   mesh_grid_t *temp_point = this->mesh_grid_pt2;
   this->mesh_grid_pt2= this->mesh_grid_pt1;
   this->mesh_grid_pt1 = temp_point; //échange des deux pointeurs
   
   this->thrust_saver();
}

void Diff_Eq_Solver::calc_iteration_PG_cyl() {
        mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
        mesh_grid_t &mesh_grid_2 = (*(this->mesh_grid_pt2));
        
        register int i,j;
        for (i = 1 ; i < this->arglist_pt->x_size-1; i++) {
                for (j = 1 ; j < this->arglist_pt->y_size-1; j++) {
      
                        this->update_vol_mass_cyl(i,j);
        
                        this->update_speed_r_cyl(i,j);
        
                        this->update_speed_z_cyl(i,j);
        
                        this->update_temp_PG_cyl(i,j);
        
                        this->update_pres_PG(i,j);
        
                }
        }
        for (i = 1 ; i < this->arglist_pt->x_size-1; i++) {
                copy_case(i,0,i,1);
        }
        for (j = 1 ; j < this->arglist_pt->y_size-1; j++) {
                copy_case(0,j,1,j);
                copy_case(this->arglist_pt->x_size-1,j,this->arglist_pt->x_size-2,j);
        }

        mesh_grid_t *temp_point = this->mesh_grid_pt2;
        this->mesh_grid_pt2 = this->mesh_grid_pt1;
        this->mesh_grid_pt1 = temp_point; //échange des deux pointeurs
}

void Diff_Eq_Solver::calc_iteration_VDW_cart() {
        mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
        mesh_grid_t &mesh_grid_2 = (*(this->mesh_grid_pt2));
        register int i,j;
        for (i = 1 ; i < this->arglist_pt->x_size-1; i++) {
                for (j = 1 ; j < this->arglist_pt->y_size-1; j++) {
      
                        this->update_vol_mass(i,j);
        
                        this->update_speed_x(i,j);
        
        this->update_speed_y(i,j);
        
        this->update_temp_VDW(i,j);
        
        this->update_pres_VDW(i,j);
        
      }
   }
   for (i = 1 ; i < this->arglist_pt->x_size-1; i++) {
           copy_case(i,0,i,1);
   }
   for (j = 1 ; j < this->arglist_pt->y_size-1; j++) {
           copy_case(0,j,1,j);
           copy_case(this->arglist_pt->x_size-1,j,this->arglist_pt->x_size-2,j);
   }
   
   mesh_grid_t *temp_point = this->mesh_grid_pt2;
   this->mesh_grid_pt2 = this->mesh_grid_pt1;
   this->mesh_grid_pt1 = temp_point; //échange des deux pointeurs
}

void Diff_Eq_Solver::calc_iteration_VDW_cyl() {
        mesh_grid_t &mesh_grid_1 = (*(this->mesh_grid_pt1));
        mesh_grid_t &mesh_grid_2 = (*(this->mesh_grid_pt2));
        
        register int i,j;
        for (i = 1 ; i < this->arglist_pt->x_size-1; i++) {
                for (j = 1 ; j < this->arglist_pt->y_size-1; j++) {
      
                        this->update_vol_mass_cyl(i,j);
      

                        this->update_speed_r_cyl(i,j);
        
                        this->update_speed_z_cyl(i,j);
        
                        this->update_temp_VDW_cyl(i,j);
        
                        this->update_pres_VDW(i,j);
        
                }
        }
        for (i = 1 ; i < this->arglist_pt->x_size-1; i++) {
                copy_case(i,0,i,1);
        }
        for (j = 1 ; j < this->arglist_pt->y_size-1; j++) {
                copy_case(0,j,1,j);
                copy_case(this->arglist_pt->x_size-1,j,this->arglist_pt->x_size-2,j);
        }

        mesh_grid_t *temp_point = this->mesh_grid_pt2;
        this->mesh_grid_pt2 = this->mesh_grid_pt1;
        this->mesh_grid_pt1 = temp_point; //échange des deux pointeurs
}


void Diff_Eq_Solver::solve_PG_cart()
{
        register int i;
        for (i=0; i<this->arglist_pt->iter_number_solver; i++) { this->calc_iteration_PG_cart(); }
}

void Diff_Eq_Solver::solve_VDW_cart()
{
        register int i;
        for (i=0; i<this->arglist_pt->iter_number_solver; i++) { this->calc_iteration_VDW_cart(); }
}

void Diff_Eq_Solver::solve_PG_cyl()
{
        register int i;
        for (i=0; i<this->arglist_pt->iter_number_solver; i++) { this->calc_iteration_PG_cyl(); }
}

void Diff_Eq_Solver::solve_VDW_cyl()
{
        register int i;
        for (i=0; i<this->arglist_pt->iter_number_solver; i++) { this->calc_iteration_VDW_cyl(); }
}


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

                default:
                        throw "DES : Bad solving algorithm";

        }
}
