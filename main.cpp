/* This is main file 
 *
 * Coded by J.Auffiger; A.Valade
 *
*/


// ======== includes ========

#include "nozzle_profiler.h"
#include "data_mapper.h"
#include "usr_interface.h"
#include "arg_interpreter.h"

#include <string.h>
// ======= main code ========



int main(int argc, char *argv[])
{
       // chargement de l'interface utilisateur et de l'interpreteur d'arguments
       Usr_Interface *UI = new Usr_Interface();
       UI->cout_str("");
       UI->cout_str("------ Program launched ------");
       UI->cout_str("Loading Arguments Interpreter...");
       Arg_Interpreter *AI = new Arg_Interpreter(UI);
       
       string argfile_name = ""; //moche mais nécessaire pour le DM
       
       // gestion des arguments
              
       // cas de l'absence d'argument
       if (argc == 1) { 
               UI->cout_err("Missing arguments: type -h for help");
               return -1;
       }

       // cas de l'aide
       if (strcmp(argv[1],"-h") == 0){
               UI->cout_str_no_space("You asked for some help");
               UI->cout_str_no_space("Two options available:");
               UI->cout_str_no_space("1/ --run-config argfile.arg");
               UI->cout_str_no_space(">>> Run the given argfile");
               UI->cout_str_no_space("2/ --new-config argfile.arg");
               UI->cout_str_no_space(">>> Create a new configuration file interactively");
               return -1;
       }
       for (int i=1; i<argc; i++){

               // creation de fichiers à partir de la console; ne lance pas automatiquement le programme
               if (strcmp(argv[i],"--new-config") == 0){
                       if (i+1 < argc){
                               try{
                                       argfile_name = argv[i+1];
                                       AI->create_argfile_from_cons( argfile_name );
                                       UI->cout_str("Configuration file created -- Now exiting with no error");
                                       return 1;
                               }
                               catch (const char * err){
                                       UI->cout_err(err);
                                       return -1;
                               }
                       }
                       else {
                               UI->cout_err("Missing configuration name to create: type -h for help");
                               return -1;
                       }
               }

               // lancement du programme à partir des fichiers de conf déja écrits
               if (strcmp(argv[i],"--run-config") == 0){
                       if (i+1 < argc){
                               try{
                                       argfile_name = argv[i+1];
                                       AI->fill_arglist_from_argfile( argfile_name );
                                       break;
                               }
                               catch (const char * err){
                                       UI->cout_err(err);
                                       return -1;
                               }
                       }
                       else {
                               UI->cout_err("Missing configuration name to load: type -h for help");
                               return -1;
                       }
               }

               else {
                       UI->cout_err("Invalid argument(s): type -h for help");
                       return -1;
               }

       }
       
       UI->cout_str("Loading Data Mapper...");
       Data_Mapper *DM = new Data_Mapper(UI, AI->get_arglist_pt(), argfile_name);

       UI->cout_str("Loading Nozzle Profiler...");
       Nozzle_Profiler *NP = new Nozzle_Profiler(UI, DM, AI->get_arglist_pt());
       
       try {
               NP->profile();
       }
       catch (const char * err){
               UI->cout_err(err);
               return -1;
       }

       UI->cout_str("------ Exiting program succesfully ------");
       
      // delete AI, UI, NP;       
       
       return 1;
}
