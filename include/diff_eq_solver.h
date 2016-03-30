/* This is differential equation solver header
 *
 * This module is the one which makes the experimentation 
 *
 * The object defined below is constructed in nozzle_profiler.h and is called through its adress elsewhere
 *
*/

// ======== include guards ========
#include <vector>
#ifndef __DIFF_EQ_SOLVER_IS_INCLUDED__
#define __DIFF_EQ_SOLVER_IS_INCLUDED__

using namespace std;

// ========= obj dependencies ========
class Usr_Interface;
struct arglist_struct;
struct mesh_struct;
typedef vector<vector<mesh_struct>> mesh_grid_t;
typedef double data_t; 
// ======= includes ========


// ======== interface ========

class Diff_Eq_Solver
{
        public:
                Diff_Eq_Solver(Usr_Interface *UI,arglist_struct *arglist_pt,mesh_grid_t *mesh_grid_pt1,mesh_grid_t *mesh_grid_pt2); //constructeur

                void solve(); // méthode principale avec un switch
                
                void solve_PG_cart(); //Comme la suite mais avec les itérations temporelles !
                void solve_VDW_cart();
                void solve_PG_cyl();
                void solve_VDW_cyl();

                void calc_iteration_PG_cart(); //méthode principale gaz parfait
                void calc_iteration_VDW_cart(); //méthode principale Van der Waals
                void calc_iteration_PG_cyl(); //méthode principale gaz parfait, yoloswag en cylindriques
                void calc_iteration_VDW_cyl();//méthode principale Van der Waals en cylindriques
                void copy_case(int i, int j, int k, int l); //copie de la struct mesh_grid pour les ocnditions au bord
                
                void update_vol_mass(int i,int j); //mise à jour de la masse volumique
                void update_speed_x(int i, int j); //mise à jour de la vitesse selon x
                void update_speed_y(int i, int j); //mise à jour de la vitesse selon y
                
                void update_vol_mass_cyl(int i, int j);//mise à jour de la masse volumique en coordonnées cylindriques
                void update_speed_r_cyl(int i, int j); //mise à jour la vitesse_r en cylindriques
                void update_speed_z_cyl(int i, int j); //mise à jour la vitesse_z en cylindriques
                
                void update_temp_PG(int i, int j); //mise à jour de la température gaz parfait
                void update_temp_PG_cyl(int i, int j);//mise à jour de la température GP en cyl
                void update_pres_PG(int i, int j); //mise à jour de la pression gaz parfait
                
                void update_temp_VDW(int i, int j);//mise à jour de la température Van der Waals
                void update_temp_VDW_cyl(int i, int j);//mise à jour de la température VDW en cyl
                void update_pres_VDW(int i, int j);//mise à jour de la pression Van der Waals
                
                //quelques fonctions pour clarifier le bazar
                data_t speed2(int i, int j, int k);//la norme au carré de la vitesse
                data_t r(int i);//définit la coordonnée radiale r en fonction de i l'ordonnée dans le tableau
                data_t diver_rhov_c(int i, int j);//la divergence en cartésiennes
                data_t deriv_x_pres(int i, int j);//première composante du gradient de pression en cartésiennes
                data_t deriv_y_pres(int i, int j);//deuxième composante du gradient de pression en cartésiennes
                data_t deriv_x_rhovx(int i, int j);//dérivée par rapport à x de rho*v_x
                data_t deriv_y_rhovx(int i, int j);//dérivée par rapport à y de rho*v_x
                data_t deriv_x_rhovy(int i, int j);//dérivée par rapport à x de rho*v_y
                data_t deriv_y_rhovy(int i, int j);//dérivée par rapport à y de rho*v_y
                data_t deriv_x_rhovz(int i, int j);//dérivée par rapport à x de rho*v_z
                data_t deriv_y_rhovz(int i, int j);//dérivée par rapport à y de rho*v_z
                data_t deriv2_x_temp(int i, int j);//dérivée seconde de la température par rapport à x
                data_t deriv2_y_temp(int i, int j);//dérivée seconde de la température par rapport à y
                data_t deriv_x_temp(int i, int j);
                data_t deriv_y_temp(int i, int j);
                data_t deriv_r_temp(int i, int j);//dérivée radiale de la température
                data_t pres_tot_PG(int i, int j);//en fait c'est rho*e+P
                data_t en_tot_PG(int i, int j);//en fait c'est e l'énergie massique
                data_t pres_tot_VDW(int i, int j);//en fait c'est encore rho*e+P
                data_t en_tot_VDW(int i, int j);//en fait c'est l'énergie massique
                data_t deriv_x_prestotPGvx(int i, int j);//tout est dans le nom
                data_t deriv_y_prestotPGvy(int i, int j);//tout est dans le nom
                data_t deriv_r_prestotPGrvr(int i, int j);//tout est dans le nom
                data_t deriv_x_prestotVDWvx(int i, int j);//pareil pour VDW
                data_t deriv_y_prestotVDWvy(int i, int j);//pareil pour VDW
                data_t deriv_r_prestotVDWrvr(int i, int j);//pareil pour VDW
                data_t deriv_r_rrhovr(int i, int j);//tout est dans le nom
                data_t deriv_x_vx(int i, int j);
                data_t deriv_y_vx(int i, int j);
                data_t deriv_x_vy(int i, int j);
                data_t deriv_y_vy(int i, int j);
                data_t mol_stress_xy(int i, int j);
                data_t mol_stress_xx(int i, int j);
                data_t mol_stress_yy(int i, int j);
                data_t strain_xy(int i, int j);
                data_t strain_xx(int i, int j);
                data_t strain_yy(int i, int j);
                data_t deriv_x_tauxx(int i, int j);
                data_t deriv_x_tauyx(int i, int j);
                data_t deriv_y_tauyy(int i, int j);
                data_t deriv_y_tauxy(int i, int j);
                data_t diver_vtau(int i, int j);
                data_t vtaux(int i, int j);
                data_t vtauy(int i, int j);
                data_t deriv_x_vtaux(int i, int j);
                data_t deriv_y_vtauy(int i, int j);
                bool is_in(int i, int j);
        
        private:
                Usr_Interface *UI; //pointeur vers l'object d'interface utilisateur
                arglist_struct *arglist_pt; //pointeur vers la arglist
                mesh_grid_t *mesh_grid_pt1; //pointeur vers le tableau tuyère
                mesh_grid_t *mesh_grid_pt2;
};


#endif





