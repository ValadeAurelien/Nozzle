/* This is arguments interpreted implementation linked to arg_interpreter.h
 *
 *
 *
 *
*/


// ======== includes ======== 
#include "arg_interpreter.h"
#include "usr_interface.h"

#include <fstream>
#include <iostream>

#define READ file >> arg_name >> equal >> arg_val; this->UI->cout_str_no_space("--> " + arg_name + " = " + arg_val);
// ======== implementation ========


Arg_Interpreter::Arg_Interpreter(Usr_Interface *UI)
{
        this->UI = UI;
}

arglist_struct* Arg_Interpreter::get_arglist_pt()
{
        return &(this->arglist);
}

void Arg_Interpreter::create_argfile_from_cons(string argfile_name) {

        string argfile_path = this->PATH + argfile_name;
        this->UI->cout_str("AI-> creating argfile " + argfile_path);

        ofstream file;
        file.open(argfile_path,ios::out); //création du fichier txt arg_file

        //caractéristiques du pattern
        (this->UI)->cout_str("Enter the height of the pattern (int)");
        file << "x_size = " << (this->UI)->cin_str() << endl;

        (this->UI)->cout_str("Enter the width of the pattern (int)");
        file << "y_size = " << (this->UI)->cin_str() << endl;

        (this->UI)->cout_str("Enter the initial time step (float)");
        file << "time_step = " << (this->UI)->cin_str() << endl;

        (this->UI)->cout_str("Enter the space step (float)");
        file << "space_step = " << (this->UI)->cin_str() << endl;

        (this->UI)->cout_str("Enter the CFL condition (float)");
        file << "CFL_cond = " << (this->UI)->cin_str() << endl;

        (this->UI)->cout_str("Enter the number of iterations for the eq_diff solver (int)");
        file << "iter_number_solver = " << (this->UI)->cin_str() << endl;

        (this->UI)->cout_str("Enter the number of iterations for the profiler (int)");
        file << "iter_number_profiler = " << (this->UI)->cin_str() << endl;

        (this->UI)->cout_str("Enter the number of threads allocated to the program (int)");
        file << "nb_of_threads = " << (this->UI)->cin_str() << endl;

        //Paramètres de fitting, et conditions initiales


        this->UI->cout_str("Choose how many points you want to usei for the initial nozzle profile (int)");
        string answer = this->UI->cin_str();
        file << "->nb_pts = " << answer << endl;
                
        this->UI->cout_str("Give the abscisses of your points (int)");
        for (int i=0; i < stoi(answer); i++) {
                file << "->abs = " << this->UI->cin_str() << endl;
        }

        this->UI->cout_str("Give the ordinates of your points (int)");
        for (int i=0; i < stoi(answer); i++) {
                file << "->ord = " << this->UI->cin_str() << endl;
        }

        (this->UI)->cout_str("Choose the type of initial condition (initial_gradient=0, evoluting_chamber_conditions=1");
        answer = (this->UI)->cin_str();
        file << "init_cond_type = " << answer << endl;

        if (answer == "1"){
            (this->UI)->cout_str("Enter the number of iterations to reach the final chamber conditions (int)");
            file << "iter_number_evol_chamber = " << (this->UI)->cin_str() << endl;

        }
        else {file << "iter_number_evol_chamber = " << 0 << endl;}
        
        (this->UI)->cout_str("Choose the differential equation solver algorithm (PG_cart=0,VDW_cart=1,PG_cyl=2,VDW_cyl=3,PG_cart_turb=4)");
        answer = (this->UI)->cin_str();
        file << "differential_equation_solver_algo = " << answer << endl;
        
        // Constantes numériques pour l'eq diff

        string diff_eq_solver_algo = answer; //degueu mais nécéssaire pour plus tard --> constante pour les gaz et turbulences

        if ( answer == "1" or answer == "3" ) {
                (this->UI)->cout_str("Enter the Van der Waals coefficient 'a' (float)");
                file << "VDW_a_coeff = " << (this->UI)->cin_str() << endl;

                (this->UI)->cout_str("Enter the Van der Waals coefficient 'b' (float)");
                file << "VDW_b_coeff = " << (this->UI)->cin_str() << endl;
        }
        else {
                file << "VDW_a_coeff = 0" << endl;
                file << "VDW_b_coeff = 0" << endl;
        }

        (this->UI)->cout_str("Activate thermal conduction ? (bool: true/false)");
        answer = (this->UI)->cin_str();
        file << "thermal_conduction = " << answer << endl;

        if (answer == "true") {
        (this->UI)->cout_str("Enter the thermal conduction value (float)");
        file << "lambda = " << (this->UI)->cin_str() << endl;
        }
        else file << "lambda = 0" << endl;

        (this->UI)->cout_str("Enter the mole mass of the gaz (float)");
        file << "mol_mass = " << (this->UI)->cin_str() << endl;

        (this->UI)->cout_str("Enter the coefficient of dynamic viscosity (float)");
        file << "dyn_visc = " << (this->UI)->cin_str() << endl;

        //caractéristiques numériques initiales
        
        (this->UI)->cout_str("Enter the initial pressure of the chamber (float)");
        file << "init_chamber_pressure = " << (this->UI)->cin_str() << endl;

        (this->UI)->cout_str("Enter the initial temperature of the chamber (float)");
        file << "init_chamber_temperature = " << (this->UI)->cin_str() << endl;

        (this->UI)->cout_str("Enter the initial speed of the chamber (float)");
        file << "init_chamber_speed = " << (this->UI)->cin_str() << endl;

        if ( diff_eq_solver_algo == "4") {
                (this->UI)->cout_str("Enter the value of k (turbulence energy) in the chamber (float)");
                file << "init_chamber_turb_en = " << (this->UI)->cin_str() << endl;
                
                (this->UI)->cout_str("Enter the value of epsilon (turbulence dissipation) in the chamber (float)");
                file << "init_chamber_turb_dis = " << (this->UI)->cin_str() << endl;
        }
        else {
                file << "init_chamber_turb_en = 0" << endl;
                file << "init_chamber_turb_dis = 0" << endl;
        }

        (this->UI)->cout_str("Enter the initial pressure of the atmostphere (float)");
        file << "init_atmosphere_pressure = " << (this->UI)->cin_str() << endl;

        (this->UI)->cout_str("Enter the initial temperature of the atmostphere (float)");
        file << "init_atmosphere_temperature = " << (this->UI)->cin_str() << endl;

        (this->UI)->cout_str("Enter the initial speed of the atmosphere (float)");
        file << "init_atmosphere_speed = " << (this->UI)->cin_str() << endl;

        if ( diff_eq_solver_algo == "4") {
                (this->UI)->cout_str("Enter the value of k (turbulence energy) in the atmosphere (float)");
                file << "init_atmosphere_turb_en = " << (this->UI)->cin_str() << endl;
                
                (this->UI)->cout_str("Enter the value of epsilon (turbulence dissipation) in the atmosphere (float)");
                file << "init_atmosphere_turb_dis = " << (this->UI)->cin_str() << endl;
        }
        else {
                file << "init_atmosphere_turb_en = 0" << endl;
                file << "init_atmosphere_turb_dis = 0" << endl;
        }

        file.close();
}

void Arg_Interpreter::fill_arglist_from_argfile(string argfile_name) {
  

        string argfile_path = this->PATH + argfile_name;

        this->UI->cout_str("AI-> creating arglist from argfile " + argfile_path + " :" );
        
        ifstream file(argfile_path,ios::in);

        string equal, arg_name, arg_val;

        // caractérisation du patern
        READ 
        try {this->arglist.x_size = stoi(arg_val);}
        catch (...) {throw "AI: Invalid x size";}
    
        READ
        try {this->arglist.y_size = stoi(arg_val);}
        catch (...) {throw "AI: Invalid y size";}

        READ
        try {this->arglist.init_time_step = stod(arg_val);}
        catch (...) {throw "AI: Invalid initial time step";}

        READ
        try {this->arglist.space_step = stod(arg_val);}
        catch (...) {throw "AI: Invalid space space";}

        READ
        try {this->arglist.CFL_cond = stod(arg_val);}
        catch (...) {throw "AI: Invalid CFL condition";}

        READ
        try {this->arglist.iter_number_solver = stoi(arg_val);}
        catch (...) {throw "AI: Invalid number of iteration on differential equation solver";}

        READ
        try {this->arglist.iter_number_profiler = stoi(arg_val);}
        catch (...) {throw "AI: Invalid number of iteration on nozzle profiler";}

        READ
        try {this->arglist.nb_of_threads = stoi(arg_val);}
        catch (...) {throw "AI: Invalid number of threads";}

        // Paramètre de fitting, choix algo et arguments pour les conditions initiales

        
        READ
        try {this->arglist.nozzle_fitting_init_arg.nb_pts = stoi(arg_val);}
        catch (...) {throw "AI: Invalid number of points for the initial profile";}
        
        for (int i=0; i<this->arglist.nozzle_fitting_init_arg.nb_pts; i++){
                READ
                try {this->arglist.nozzle_fitting_init_arg.abscisses.push_back( stoi(arg_val) );}
                catch (...) {throw "AI: Bad abcisse value";}
        }

        for (int i=0; i<this->arglist.nozzle_fitting_init_arg.nb_pts; i++){
                READ
                try {this->arglist.nozzle_fitting_init_arg.ordinates.push_back( stoi(arg_val) );}
                catch (...) {throw "AI: Bad ordinate value";}
        }

        READ
        if (arg_val == "0") {this->arglist.init_cond_type = INIT_GRAD;}
        else {
            if (arg_val == "1") {this->arglist.init_cond_type = EVOL_CHAMBER;}
            else {throw "AI: Invalid choice of initial condition type";}
        }
    
        READ
        try {this->arglist.iter_number_evol_chamber = stoi(arg_val);}
        catch (...) {throw "AI: Invalid a number of iteration to reach the final chamber conditions";}
        if (this->arglist.iter_number_evol_chamber + 2 >= this->arglist.iter_number_solver) {
            throw "AI: Number of iteration to reach the final chamber conditions must be strictly inferior to the number of iteration of the solver minus o,e";
        }


        READ 
        if (arg_val == "0") {this->arglist.diff_eq_solver_algo = PG_cart;}
        else {
                if (arg_val == "1") {this->arglist.diff_eq_solver_algo = VDW_cart;}
                else {
                        if (arg_val == "2") {this->arglist.diff_eq_solver_algo = PG_cyl;}
                        else {
                                if (arg_val == "3") {this->arglist.diff_eq_solver_algo = VDW_cyl;}
                                        if (arg_val == "4") {this->arglist.diff_eq_solver_algo = PG_cart_turb;}
                                        else {throw "AI: Invalid choice of differential solver algorithm";}
                        }
                }
        }

        // Constantes numériques pour l'eq diff

        READ
        try {this->arglist.VDW_a_coef = stof(arg_val);}
        catch (...) {throw "AI: Invalid a coefficient of Van Der Waals";}

        READ
        try {this->arglist.VDW_b_coef = stof(arg_val);}
        catch (...) {throw "AI: Invalid b coefficient of Van Der Waals";}
    
        READ
        if (arg_val == "false")  {this->arglist.thermal_conduction = false;}
        else {
                if (arg_val == "true")  {this->arglist.thermal_conduction = true;}
                else {throw "AI: Invalid choice of thermal conduction";}
        }

        READ
        try {this->arglist.lambda = stof(arg_val);}
        catch (...) {throw "AI: Invalid lambda";}

        //caractéristiques numériques initiales

        READ
        try {this->arglist.mol_mass = stof(arg_val);}
        catch (...) {throw "AI: Invalid molar mass";}

        READ
        try {this->arglist.dyn_visc = stof(arg_val);}
        catch (...) {throw "AI: Invalid dynamic viscosity";}

        READ
        try {this->arglist.init_cond.chamber_pressure = stof(arg_val);}
        catch (...) {throw "AI: Invalid initial chamber pressure";}

        READ
        try {this->arglist.init_cond.chamber_temp = stof(arg_val);}
        catch (...) {throw "AI: Invalid initial chamber temperature";}

        READ
        try {this->arglist.init_cond.chamber_speed = stof(arg_val);}
        catch (...) {throw "AI: Invalid initial chamber gaz speed";}

        READ
        try {this->arglist.init_cond.chamber_turb_en = stof(arg_val);}
        catch (...) {throw "AI: Invalid value of initial chamber turbulence energy";}
        
        READ
        try {this->arglist.init_cond.chamber_turb_dis = stof(arg_val);}
        catch (...) {throw "AI: Invalid value of initial chamber turbulence dissipation";}
        
        READ
        try {this->arglist.init_cond.atmosphere_pressure = stof(arg_val);}
        catch (...) {throw "AI: Invalid initial atmosphere pressure";}

        READ
        try {this->arglist.init_cond.atmosphere_temp = stof(arg_val);}
        catch (...) {throw "AI: Invalid initial atmosphere temperature";}

        READ
        try {this->arglist.init_cond.atmosphere_speed = stof(arg_val);}
        catch (...) {throw "AI: Invalid initial atmosphere gaz speed";}

        READ
        try {this->arglist.init_cond.atmosphere_turb_en = stof(arg_val);}
        catch (...) {throw "AI: Invalid value of initial atmosphere turbulence energy";}
        
        READ
        try {this->arglist.init_cond.atmosphere_turb_dis = stof(arg_val);}
        catch (...) {throw "AI: Invalid value of initial atmosphere turbulence dissipation";}

        this->UI->new_line();

        file.close();
}
