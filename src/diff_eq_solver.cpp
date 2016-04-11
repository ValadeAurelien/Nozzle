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
#include <cmath>
#include <vector>

#define R 8.314
#define sig_k 1
#define sig_e 1.3
#define c_e_1 1.45
#define c_e_2 1.92
#define c_mu 0.09
#define sig_t 1.3

#define mesh_grid_1 (*(this->mesh_grid_pt1))
#define mesh_grid_2 (*(this->mesh_grid_pt2))

// ======== implementation ========

Diff_Eq_Solver::Diff_Eq_Solver(Usr_Interface *UI,Data_Mapper *DM, arglist_struct *arglist_pt, mesh_grid_t *mesh_grid_pt1, mesh_grid_t *mesh_grid_pt2) {
        this->UI = UI;
        this->DM = DM;
        this->arglist_pt = arglist_pt;
        this->mesh_grid_pt1 = mesh_grid_pt1;
        this->mesh_grid_pt2 = mesh_grid_pt2;
        vector<data_t> thrust(this->arglist_pt->iter_number_solver);
        this->thrust = thrust;
}

data_t Diff_Eq_Solver::speed2(int i, int j, int k) {
        
        if (k==1) {
                return(pow(mesh_grid_1[i][j].speed[0],2) + pow(mesh_grid_1[i][j].speed[1],2));
        }
        else {
                return(pow(mesh_grid_2[i][j].speed[0],2) + pow(mesh_grid_2[i][j].speed[1],2));
        }
}

data_t Diff_Eq_Solver::pres_tot_PG(int i, int j) {
        
        return(
        mesh_grid_1[i][j].vol_mass*en_tot_PG(i,j)+mesh_grid_1[i][j].pressure
        );
}

data_t Diff_Eq_Solver::pres_tot_VDW(int i, int j) {
        
        return(
        mesh_grid_1[i][j].vol_mass*en_tot_VDW(i,j)+mesh_grid_1[i][j].pressure
        );
}

data_t Diff_Eq_Solver::en_tot_PG(int i, int j) {
        
        return(
        5.0/2.0*R/this->arglist_pt->mol_mass*mesh_grid_1[i][j].temperature + 1.0/2.0*speed2(i,j,1)
        );
}

data_t Diff_Eq_Solver::en_tot_VDW(int i, int j) {
        
        return(
        5.0/2.0*R*mesh_grid_1[i][j].temperature/this->arglist_pt->mol_mass + 1.0/2.0*speed2(i,j,1) - this->arglist_pt->VDW_a_coef*mesh_grid_1[i][j].vol_mass/pow(this->arglist_pt->mol_mass,2)
        );
}

data_t Diff_Eq_Solver::diver_rhov_c(int i, int j) {
        
        return(
             deriv_x_rhovx(i,j)+deriv_y_rhovy(i,j)
        );
}

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

data_t Diff_Eq_Solver::deriv_x_rhovx(int i, int j) {
        
        return(
        (mesh_grid_1[i+1][j].vol_mass*mesh_grid_1[i+1][j].speed[0]-mesh_grid_1[i-1][j].vol_mass*mesh_grid_1[i-1][j].speed[0])
        / (2*this->arglist_pt->space_step)
        );
}

data_t Diff_Eq_Solver::deriv_y_rhovx(int i, int j) {
        
        return(
        (mesh_grid_1[i][j+1].vol_mass*mesh_grid_1[i][j+1].speed[0]-mesh_grid_1[i][j-1].vol_mass*mesh_grid_1[i][j-1].speed[0])
        / (this->arglist_pt->space_step*2)
        );
}

data_t Diff_Eq_Solver::deriv_x_rhovy(int i, int j) {
        
        return(
        (mesh_grid_1[i+1][j].vol_mass*mesh_grid_1[i+1][j].speed[1]-mesh_grid_1[i-1][j].vol_mass*mesh_grid_1[i-1][j].speed[1])
        / (this->arglist_pt->space_step*2)
        );
}

data_t Diff_Eq_Solver::deriv_y_rhovy(int i, int j) {
        
        return(
        (mesh_grid_1[i][j+1].vol_mass*mesh_grid_1[i][j+1].speed[1]-mesh_grid_1[i][j-1].vol_mass*mesh_grid_1[i][j-1].speed[1])
        / (this->arglist_pt->space_step*2)
        );
}

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

data_t Diff_Eq_Solver::deriv_r_prestotPGrvr(int i, int j) {
        
        return(
        -(r(i+1)*mesh_grid_1[i+1][j].speed[0]*pres_tot_PG(i+1,j)-r(i-1)*mesh_grid_1[i-1][j].speed[0]*pres_tot_PG(i-1,j))
        / (this->arglist_pt->space_step*2)
        );
}

data_t Diff_Eq_Solver::deriv_r_temp(int i, int j) {
        
        return(
        (1-mesh_grid_1[i+1][j].is_wall)*(1-mesh_grid_1[i-1][j].is_wall)*(mesh_grid_1[i-1][j].temperature-mesh_grid_1[i+1][j].temperature)
        / (this->arglist_pt->space_step*2)
        );
}

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

data_t Diff_Eq_Solver::deriv_r_prestotVDWrvr(int i, int j) {
        
        return(
        -(r(i+1)*mesh_grid_1[i+1][j].speed[0]*pres_tot_VDW(i+1,j)-r(i-1)*mesh_grid_1[i-1][j].speed[0]*pres_tot_VDW(i-1,j))
        / (this->arglist_pt->space_step*2)
        );
}

data_t Diff_Eq_Solver::r(int i) {
        return(this->arglist_pt->space_step*(this->arglist_pt->x_size-i));
}

data_t Diff_Eq_Solver::deriv_r_rrhovr(int i, int j) {
        
        return(
        (r(i-1)*mesh_grid_1[i-1][j].vol_mass*mesh_grid_1[i-1][j].speed[0]-r(i+1)*mesh_grid_1[i+1][j].vol_mass*mesh_grid_1[i+1][j].speed[0])
        / (this->arglist_pt->space_step*2)
        );
}

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
	return(
	deriv_x_vtaux(i,j)+deriv_y_vtauy(i,j)
	);
}

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

data_t Diff_Eq_Solver::mu_t(int i, int j) {
	
	return( c_mu*mesh_grid_1[i][j].vol_mass*pow(mesh_grid_1[i][j].turb_en,2)/mesh_grid_1[i][j].turb_dis);
}

data_t Diff_Eq_Solver::lambda_t(int i, int j) {
	return(
	5.0/2.0*R/this->arglist_pt->mol_mass*mu_t(i,j)/sig_t
	);
}

data_t Diff_Eq_Solver::heat_flux_x_turb(int i, int j) {
	return((this->arglist_pt->lambda+lambda_t(i,j))*deriv_x_temp(i,j));
}

data_t Diff_Eq_Solver::heat_flux_y_turb(int i, int j) {
	return((this->arglist_pt->lambda+lambda_t(i,j))*deriv_y_temp(i,j));
}

data_t Diff_Eq_Solver::diver_heat_flux_turb(int i, int j) {
	return(deriv_x_heat_flux_x_turb(i,j)+deriv_y_heat_flux_y_turb(i,j));
}

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
		(heat_flux_y_turb(i+1,j)-heat_flux_y_turb(i-1,j))
		/ (this->arglist_pt->space_step*2)
		);
	}
	else {
	return(0);
	}
}

data_t Diff_Eq_Solver::diver_vtau_turb(int i, int j) {
	return(
	deriv_x_vtaux_turb(i,j)+deriv_y_vtauy_turb(i,j)
	);
}

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

data_t Diff_Eq_Solver::tot_stress_yy(int i, int j) {
	return(mol_stress_yy(i,j)+turb_stress_yy(i,j));
}

data_t Diff_Eq_Solver::turb_stress_yy(int i, int j) {
	return(2.0*mu_t(i,j)*(2.0/3.0*strain_yy(i,j)-1.0/3.0*strain_xx(i,j))-2.0/3.0*mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].turb_en);
}

data_t Diff_Eq_Solver::tot_stress_xx(int i, int j) {
	return(mol_stress_xx(i,j)+turb_stress_xx(i,j));
}

data_t Diff_Eq_Solver::turb_stress_xx(int i, int j) {
	
	return(2.0*mu_t(i,j)*(2.0/3.0*strain_xx(i,j)-1.0/3.0*strain_yy(i,j))-2.0/3.0*mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].turb_en);
}

data_t Diff_Eq_Solver::tot_stress_xy(int i, int j) {
	return(mol_stress_xy(i,j)+turb_stress_xy(i,j));
}

data_t Diff_Eq_Solver::turb_stress_xy(int i, int j) {
	return(2.0*mu_t(i,j)*strain_xy(i,j));
}

data_t Diff_Eq_Solver::diver_rhovk(int i, int j) {
	return(deriv_x_rhovxk(i,j)+deriv_y_rhovyk(i,j));
}

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

data_t Diff_Eq_Solver::diver_mudk(int i, int j) {
	return(deriv_x_mudxk(i,j)+deriv_y_mudyk(i,j));
}

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

data_t Diff_Eq_Solver::deriv_x_k(int i, int j) {
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

data_t Diff_Eq_Solver::diver_rhovepsilon(int i, int j) {
	return(deriv_x_rhovxepsilon(i,j)+deriv_y_rhovyepsilon(i,j));
}

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

data_t Diff_Eq_Solver::diver_mudepsilon(int i, int j) {
	return(deriv_x_mudxepsilon(i,j)+deriv_y_mudyepsilon(i,j));
}

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

data_t Diff_Eq_Solver::deriv_x_epsilon(int i, int j) {
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

bool Diff_Eq_Solver::is_in(int i, int j) {
	if (i>=0 && i<this->arglist_pt->x_size && j>=0 && j<this->arglist_pt->y_size) {
		return(true);
	}
	else {
		return(false);
	}
}

void Diff_Eq_Solver::update_vol_mass(int i, int j) { //mise à jour de la masse volumique
        
        if (not (mesh_grid_1[i][j].is_wall)) {
                mesh_grid_2[i][j].vol_mass = mesh_grid_1[i][j].vol_mass - diver_rhov_c(i,j)*this->arglist_pt->time_step;
        }
}

void Diff_Eq_Solver::update_vol_mass_cyl(int i, int j) {
        
        if (not (mesh_grid_1[i][j].is_wall)){
                mesh_grid_2[i][j].vol_mass = mesh_grid_1[i][j].vol_mass-this->arglist_pt->time_step*(1.0/r(i)*deriv_r_rrhovr(i,j)+deriv_y_rhovy(i,j));
        }
}

void Diff_Eq_Solver::update_speed_x(int i, int j) {
	
	if (not(mesh_grid_1[i][j].is_wall)) {
	mesh_grid_2[i][j].speed[0] = 1.0/mesh_grid_2[i][j].vol_mass*(mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].speed[0]-this->arglist_pt->time_step*(mesh_grid_1[i][j].speed[0]*diver_rhov_c(i,j)+mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].speed[0]*deriv_x_vx(i,j)+mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].speed[1]*deriv_y_vx(i,j)+deriv_x_pres(i,j)-deriv_y_tauxy(i,j)-deriv_x_tauxx(i,j) ) );
	}
}

void Diff_Eq_Solver::update_speed_y(int i, int j) {

	if (not(mesh_grid_1[i][j].is_wall)) {
	mesh_grid_2[i][j].speed[1] = 1.0/mesh_grid_2[i][j].vol_mass*(mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].speed[1]-this->arglist_pt->time_step*(mesh_grid_1[i][j].speed[1]*diver_rhov_c(i,j)+mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].speed[0]*deriv_x_vy(i,j)+mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].speed[1]*deriv_y_vy(i,j)+deriv_y_pres(i,j)-deriv_x_tauyx(i,j)-deriv_y_tauyy(i,j)  )  );
	}
}

void Diff_Eq_Solver::update_speed_r_cyl(int i, int j) {
        
        if (not (mesh_grid_1[i][j].is_wall) ){
                mesh_grid_2[i][j].speed[0] = 1.0/mesh_grid_2[i][j].vol_mass*(mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].speed[0]-this->arglist_pt->time_step*(-deriv_x_pres(i,j)-mesh_grid_1[i][j].speed[0]*deriv_x_rhovx(i,j)+mesh_grid_1[i][j].speed[1]*deriv_y_rhovy(i,j)));
        }
}

void Diff_Eq_Solver::update_speed_z_cyl(int i, int j) {
        
        if (not (mesh_grid_1[i][j].is_wall)) {
                mesh_grid_2[i][j].speed[1] = 1.0/mesh_grid_2[i][j].vol_mass*(mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].speed[1]-this->arglist_pt->time_step*(deriv_y_pres(i,j) - mesh_grid_1[i][j].speed[0]*deriv_x_rhovx(i,j)+mesh_grid_1[i][j].speed[1]*deriv_y_rhovy(i,j) )  );
        }
}

void Diff_Eq_Solver::update_speed_x_turb(int i, int j) {

	if (not(mesh_grid_1[i][j].is_wall)) {
	mesh_grid_2[i][j].speed[0] = 1.0/mesh_grid_2[i][j].vol_mass*(mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].speed[0]-this->arglist_pt->time_step*(mesh_grid_1[i][j].speed[0]*diver_rhov_c(i,j)+mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].speed[0]*deriv_x_vx(i,j)+mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].speed[1]*deriv_y_vx(i,j)+deriv_x_pres(i,j)-deriv_y_tauxy_turb(i,j)-deriv_x_tauxx_turb(i,j) ) );
	}
}

void Diff_Eq_Solver::update_speed_y_turb(int i, int j) {

	if (not(mesh_grid_1[i][j].is_wall)) {
	mesh_grid_2[i][j].speed[1] = 1.0/mesh_grid_2[i][j].vol_mass*(mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].speed[1]-this->arglist_pt->time_step*(mesh_grid_1[i][j].speed[1]*diver_rhov_c(i,j)+mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].speed[0]*deriv_x_vy(i,j)+mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].speed[1]*deriv_y_vy(i,j)+deriv_y_pres(i,j)-deriv_x_tauyx_turb(i,j)-deriv_y_tauyy_turb(i,j)  )  );
	}
}

void Diff_Eq_Solver::update_pres_VDW(int i, int j) { //mise à jour de la pression Van der Waals

        if (not (mesh_grid_1[i][j].is_wall)) {
        mesh_grid_2[i][j].pressure = mesh_grid_2[i][j].vol_mass * R * mesh_grid_2[i][j].temperature / (this->arglist_pt->mol_mass - this->arglist_pt->VDW_b_coef * mesh_grid_2[i][j].vol_mass) - this->arglist_pt->VDW_a_coef * pow(mesh_grid_2[i][j].vol_mass / this->arglist_pt->mol_mass,2);
        }
}

void Diff_Eq_Solver::update_temp_VDW(int i, int j) { //mise à jour de la température Van der Waals

        if (not (mesh_grid_1[i][j].is_wall)) {
                mesh_grid_2[i][j].temperature = 2.0/5.0*this->arglist_pt->mol_mass/R*(this->arglist_pt->VDW_a_coef*mesh_grid_2[i][j].vol_mass/pow(this->arglist_pt->mol_mass,2)-1.0/2.0*speed2(i,j,2)+1.0/mesh_grid_2[i][j].vol_mass*(mesh_grid_1[i][j].vol_mass*en_tot_VDW(i,j)+this->arglist_pt->time_step*(this->arglist_pt->lambda*(deriv2_x_temp(i,j)+deriv2_y_temp(i,j))-deriv_x_prestotVDWvx(i,j)-deriv_y_prestotVDWvy(i,j))));
        }
}

void Diff_Eq_Solver::update_temp_VDW_cyl(int i, int j) {
        
        if (not (mesh_grid_1[i][j].is_wall) ){
                mesh_grid_2[i][j].temperature = 2.0/5.0*this->arglist_pt->mol_mass/R*(this->arglist_pt->VDW_a_coef*mesh_grid_2[i][j].vol_mass/pow(this->arglist_pt->mol_mass,2)-1.0/2.0*speed2(i,j,2)+1.0/mesh_grid_2[i][j].vol_mass*(mesh_grid_1[i][j].vol_mass*en_tot_VDW(i,j)+this->arglist_pt->time_step*(this->arglist_pt->lambda*(-deriv2_x_temp(i,j)+1.0/r(i)*deriv_r_temp(i,j)+deriv2_y_temp(i,j))-1.0/r(i)*deriv_r_prestotVDWrvr(i,j)-deriv_y_prestotVDWvy(i,j))));   
        }
}

void Diff_Eq_Solver::update_temp_PG(int i, int j) { //mise à jour de la température gaz parfait

        if (not (mesh_grid_1[i][j].is_wall) ){
                mesh_grid_2[i][j].temperature = 2.0/5.0*this->arglist_pt->mol_mass/R*(-1.0/2.0*speed2(i,j,2)+1.0/mesh_grid_2[i][j].vol_mass*(mesh_grid_1[i][j].vol_mass*en_tot_PG(i,j)+this->arglist_pt->time_step*(this->arglist_pt->lambda*(deriv2_x_temp(i,j)+deriv2_y_temp(i,j))-deriv_x_prestotPGvx(i,j)-deriv_y_prestotPGvy(i,j)+diver_vtau(i,j))));
        }
}

void Diff_Eq_Solver::update_temp_PG_cyl(int i, int j) {
        
        if (not (mesh_grid_1[i][j].is_wall)) {
                mesh_grid_2[i][j].temperature = 2.0/5.0*this->arglist_pt->mol_mass/R*(-1.0/2.0*speed2(i,j,2)+1.0/mesh_grid_2[i][j].vol_mass*(mesh_grid_1[i][j].vol_mass*en_tot_PG(i,j)+this->arglist_pt->time_step*(this->arglist_pt->lambda*(-deriv2_x_temp(i,j)+1.0/r(i)*deriv_r_temp(i,j)+deriv2_y_temp(i,j))-1.0/r(i)*deriv_r_prestotPGrvr(i,j)-deriv_y_prestotPGvy(i,j))));
        }
}

void Diff_Eq_Solver::update_temp_PG_turb(int i, int j) { //mise à jour de la température gaz parfait

        if (not (mesh_grid_1[i][j].is_wall) ){
                mesh_grid_2[i][j].temperature = 2.0/5.0*this->arglist_pt->mol_mass/R*(-1.0/2.0*speed2(i,j,2)+1.0/mesh_grid_2[i][j].vol_mass*(mesh_grid_1[i][j].vol_mass*en_tot_PG(i,j)+this->arglist_pt->time_step*(diver_heat_flux_turb(i,j)-deriv_x_prestotPGvx(i,j)-deriv_y_prestotPGvy(i,j)+diver_vtau_turb(i,j))));
        }
}

void Diff_Eq_Solver::update_pres_PG(int i, int j) { //mise à jour de la pression gaz parfait
        
        if (not (mesh_grid_1[i][j].is_wall)) {
                mesh_grid_2[i][j].pressure = mesh_grid_2[i][j].vol_mass * R * mesh_grid_2[i][j].temperature / this->arglist_pt->mol_mass;
        }
}

void Diff_Eq_Solver::update_k_PG_turb(int i, int j) {
	if (not(mesh_grid_1[i][j].is_wall)) {
		mesh_grid_2[i][j].turb_en = 1./mesh_grid_2[i][j].vol_mass*(mesh_grid_1[i][j].turb_en*mesh_grid_1[i][j].vol_mass+this->arglist_pt->time_step*(turb_stress_xx(i,j)*strain_xx(i,j)+turb_stress_yy(i,j)*strain_yy(i,j)+2*turb_stress_xy(i,j)*strain_xy(i,j)-mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].turb_dis-diver_rhovk(i,j)+diver_mudk(i,j)));
	}
}

void Diff_Eq_Solver::update_epsilon_PG_turb(int i, int j) {
	if (not(mesh_grid_1[i][j].is_wall)) {
		mesh_grid_2[i][j].turb_dis = 1./mesh_grid_2[i][j].vol_mass*(mesh_grid_1[i][j].vol_mass*mesh_grid_1[i][j].turb_dis+this->arglist_pt->time_step*(c_e_1*mesh_grid_1[i][j].turb_dis/mesh_grid_1[i][j].turb_en*(turb_stress_xx(i,j)*strain_xx(i,j)+turb_stress_yy(i,j)*strain_yy(i,j)+2*turb_stress_xy(i,j)*strain_xy(i,j)) - c_e_2*mesh_grid_1[i][j].vol_mass*pow(mesh_grid_1[i][j].turb_dis,2)/mesh_grid_1[i][j].turb_en - diver_rhovepsilon(i,j) + diver_mudepsilon(i,j)));
	}
}

void Diff_Eq_Solver::copy_case(int i, int j, int k, int l) {
        
        mesh_grid_2[i][j].vol_mass = mesh_grid_2[k][l].vol_mass;
        mesh_grid_2[i][j].temperature = mesh_grid_2[k][l].temperature;
        mesh_grid_2[i][j].pressure = mesh_grid_2[k][l].pressure;
        mesh_grid_2[i][j].speed[0] = mesh_grid_2[k][l].speed[0];
        mesh_grid_2[i][j].speed[1] = mesh_grid_2[k][l].speed[1];
}




// Une iteration totale

void Diff_Eq_Solver::exchange_mesh_grid_pts()
{
   mesh_grid_t *temp_point = this->mesh_grid_pt2;
   this->mesh_grid_pt2= this->mesh_grid_pt1;
   this->mesh_grid_pt1 = temp_point; //échange des deux pointeurs
}


void Diff_Eq_Solver::calc_iteration_PG_cart() {
        
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
   
   this->exchange_mesh_grid_pts();
}

void Diff_Eq_Solver::calc_iteration_PG_cart_turb() {
        
  register int i,j;
  for (i = 1 ; i < this->arglist_pt->x_size-1; i++) {
    for (j = 1 ; j < this->arglist_pt->y_size-1; j++) {
      
        this->update_vol_mass(i,j);
        
	this->update_speed_x_turb(i,j);
        
        this->update_speed_y_turb(i,j);
        
        this->update_temp_PG_turb(i,j);
        
        this->update_pres_PG(i,j);
        
        this->update_k_PG_turb(i,j);
        
        this->update_epsilon_PG_turb(i,j);
        
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

void Diff_Eq_Solver::calc_iteration_PG_cyl() {
        
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

        this->exchange_mesh_grid_pts();
}

void Diff_Eq_Solver::calc_iteration_VDW_cart() {
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
   
        this->exchange_mesh_grid_pts();
}

void Diff_Eq_Solver::calc_iteration_VDW_cyl() {
        
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

        this->exchange_mesh_grid_pts();
}


void Diff_Eq_Solver::solve_PG_cart()
{
	data_t thrusty = 0.;
        register int i;
        for (i=0; i<this->arglist_pt->iter_number_solver; i++) {
        	this->calc_iteration_PG_cart();
        	
        	thrusty=0;
   		for (int k = 0; k<this->arglist_pt->x_size; k++) {
      			thrusty+= -mesh_grid_1[k][0].speed[1]*pow(this->arglist_pt->space_step,2);
   		}
   		this->thrust[i] = thrusty;
        }
        this->DM->thrust_plotter(&(this->thrust));
}

void Diff_Eq_Solver::solve_PG_cart_turb() {
	data_t thrusty = 0.;
        register int i;
        for (i=0; i<this->arglist_pt->iter_number_solver; i++) {
        	this->calc_iteration_PG_cart_turb();
        	
        	thrusty=0;
   		for (int k = 0; k<this->arglist_pt->x_size; k++) {
      			thrusty+= -mesh_grid_1[k][0].speed[1]*pow(this->arglist_pt->space_step,2);
   		}
   		this->thrust[i] = thrusty;
        }
        this->DM->thrust_plotter(&(this->thrust));
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
                
                case PG_cart_turb:
                	this->solve_PG_cart_turb();
                	break;

                default:
                        throw "DES : Bad solving algorithm";

        }
}
