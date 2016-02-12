/* This is arguments files interpreter header
 * 
 * This module is the one which reads arg files and transforms thoses into c++ objects
 * It is also able to write arg files while discussing with the user in the console 
 *
 * The object defined below is constructed in main.h 
 *
*/

// ======== include guards ========

#ifndef __ARG_INTERPRETER_IS_INCLUDED__
#define __ARG_INTERPRETER_IS_INCLUDED__

// ========= obj dependencies ========

Usr_Interface *UI;

// ======= includes ========
 
#include <string>

// ======== interface ========

using namespace std;

// enums

enum gaz_type_en {GP=0, VDW=1}; //liste des gazs possibles, on assigne un entier en plus pour pouvoir des switch si on ajoute d'autres types de gaz
enum nozzle_fitting_algo_en {LAGRANGE=0}; // liste des choix de fiitting possible pour la tuyère. Utilisé dans nozzle_profiler.h

// - liste des conditions aux limites applicables
enum limit_pressure_conditions_en {};
enum limit_temperature_conditions_en {};
enum limit_vol_mass_conditions_en {};
enum limit_speed_conditions_en {NUL_SPEED_ON_BOARDERS=0};
enum init_conditions_type_en {CLASSIC=0};


// structures

struct init_conditions {
        float init_chamber_pressure;
        float init_atmosphere_pressure;
        
        float init_chamber_temp;
        float init_atmosphere_temp;

        float init_chamber_vol_mass;
        float init_atmosphere_vol_mass;

        float init_chamber_speed;
        float init_atmosphere_speed;
};

struct arg_list {
	unsigned int x_size; //le nombre de cases suivant x
	unsigned int y_size; //le nombre de cases suivant y
	
	float time_step; //le pas de temps, qui sera utilisé par eq_diff
       	float spac_step; // pas d'espace

	nozzle_fitting_algo_en nozzle_fitting_algo; //profil d'initialisation de la tuyère, qui sera utilisé par le profiler
	
	bool thermal_conduction; //vaudra true ou false selon le modèle ?
	float lambda; //la conductivité thermique (qui pourra être utile selon le modèle ?
	
	gaz_type_en gaz_type; //vaudra "GP" ou "VDW" selon le type de gaz choisi
	float VDW_a_coef; //le premier coef de Van der Waals
	float VDW_b_coef; //le deuxième coef de Van der Waals
	float PG_const; //la constante des gaz parfaits
	float mol_mass; //la masse molaire du gaz, utilisée dans l'équation d'état
	float dyn_visc; //la viscosité dynamique, au cas où on décide de l'utiliser finalement
};

// objets

class Arg_Interpreter {
        public:
		Arg_Interpreter(Usr_Interface *UI); //constructeur
		arg_list create_arglist(); //la méthode principale
                void create_argfile_from_cons();
        private:
		const string PATH = "./argfiles";
};              Usr_Interface *UI;

#endif





