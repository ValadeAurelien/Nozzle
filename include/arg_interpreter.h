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

class Usr_Interface;

// ======= includes ========
 
#include <string>
#include <vector> // pour les arguments pour l'initialissation du fitting de la tuyère
#define LIMIT_SPEED_CONDITIONS_EN {NUL_SPEED_ON_BOARDERS=0}
#define INIT_CONDITIONS_TYPE_EN {CLASSIC=0}


// ======== interface ========

using namespace std;

// enums

enum gaz_type_en {PG=0, VDW=1}; //liste des gazs possibles, on assigne un entier en plus pour pouvoir des switch si on ajoute d'autres types de gaz
//enum nozzle_fitting_algo_en {LAGRANGE=0, BRUTAL_FORCE=1, SEGMENT=2}; // liste des choix de fiitting possible pour la tuyère. Utilisé dans nozzle_profiler.h
enum init_cond_type_en {INIT_GRAD=0, EVOL_CHAMBER=1}; // gradient initial ou evolution des conditions au fond de la chambre
enum diff_eq_solver_algo_en {PG_cart=0, VDW_cart=1, PG_cyl=2, VDW_cyl=3, PG_cart_turb=4};


// structures et unions

struct init_conditions_struct {
        float chamber_pressure;
        float atmosphere_pressure;
        
        float chamber_temp;
        float atmosphere_temp;

        float chamber_speed;
        float atmosphere_speed;
        
        float chamber_turb_en;
        float atmosphere_turb_en;

        float chamber_turb_dis;
        float atmosphere_turb_dis;
};

struct nozzle_fitting_init_arg_struct {
        int nb_pts;
        vector<int> abscisses;
        vector<int> ordinates;
};

struct arglist_struct {
        unsigned int x_size; //le nombre de cases suivant x
        unsigned int y_size; //le nombre de cases suivant y
        unsigned int iter_number_solver; //le nombre d'itérations du solver d'eq_diff
        unsigned int iter_number_profiler; //le nombre d'itérations du profiler
        unsigned int nb_of_threads; //nb de threads alloués au programme

        double init_time_step; //le pas de temps, qui sera utilisé par eq_diff
        double space_step; // pas d'espace
        double CFL_cond; //rapport pas de temps/ espace et vitesse max

//        nozzle_fitting_algo_en nozzle_fitting_algo; //profil d'initialisation de la tuyère, qui sera utilisé par le profiler
        nozzle_fitting_init_arg_struct nozzle_fitting_init_arg; // arguments nécéssaires à l'initialisation de l'aglo de fitting

        init_cond_type_en init_cond_type; // la forme des conditions initiales (gradient ou evolution de la chmabre)
        unsigned int iter_number_evol_chamber;  // nb d'itérations pour mettre les conditions dans la chambre a la valeur donnée

        diff_eq_solver_algo_en diff_eq_solver_algo; //type de solveur pour le solveur d'equation différentielle

        bool thermal_conduction; //vaudra true ou false selon le modèle ?
        float lambda; //la conductivité thermique (qui pourra être utile selon le modèle ?
	
        gaz_type_en gaz_type; //vaudra "GP" ou "VDW" selon le type de gaz choisi
        float VDW_a_coef; //le premier coef de Van der Waals
        float VDW_b_coef; //le deuxième coef de Van der Waals
        float mol_mass; //la masse molaire du gaz, utilisée dans l'équation d'état
        float dyn_visc; //la viscosité dynamique, au cas où on décide de l'utiliser finalement

        init_conditions_struct init_cond; //les conditions initiales
};



// objets

class Arg_Interpreter {
        public:
		Arg_Interpreter(Usr_Interface *UI); //constructeur
		void fill_arglist_from_argfile(string argfile_name); // construction de l'objet c++ à partir d'un fichier texte formaté
                void create_argfile_from_cons(string argfile_name); // contruction d'un fichier texte structuré à partir de commandes consoles

                arglist_struct* get_arglist_pt();

        private:
		const string PATH = "./argfiles/"; // path to argfiles
                arglist_struct arglist;
                Usr_Interface *UI; //pointeur vers l'object d'interface utilisateur
};
#endif




