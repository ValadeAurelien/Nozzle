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


// ======= includes ========
 
#include <string>

// ======== interface ========

using namespace std;

struct arg_list {
	unsigned int x_size; //le nombre de cases suivant x
	unsigned int y_size; //le nombre de cases suivant y
	
	float time_step; //le pas de temps, qui sera utilisé par eq_diff
       	float spac_step; // pas d'espace

	string type_profile_init; //profil d'initialisation de la tuyère, qui sera utilisé par le profiler
	
	bool thermal_conduction; //vaudra true ou false selon le modèle ?
	float lambda; //la conductivité thermique (qui pourra être utile selon le modèle ?
	
	string type_gaz; //vaudra "GP" ou "VDW" selon le type de gaz choisi
	float VDW_a_coef; //le premier coef de Van der Waals
	float VDW_b_coef; //le deuxième coef de Van der Waals
	float PG_const; //la constante des gaz parfaits
	float mol_mass; //la masse molaire du gaz, utilisée dans l'équation d'état
	float dyn_visc; //la viscosité dynamique, au cas où on décide de l'utiliser finalement
};

class Arg_Interpreter {
        public:
		Arg_Interpreter(); //constructeur
		arg_list create_arglist(); //la méthode principale
                void create_argfile_from_cons();
        private:
		const string PATH = "./argfiles";
};

#endif





