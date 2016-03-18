
/* This is console out stream handler header
 *
 * This module is made in order control every console printing needed through out the program 
 *
 * The class defined below is created in main.h and addressed through its adress elsewhere
 *
*/

// ========= include guards ========

#ifndef __USR_INTERFACE_IS_INCLUDED__ 
#define __USR_INTERFACE_IS_INCLUDED__

// ======== obj dependencies ========


// ======== includes ========

#include <string>

// ======== interface ========

using namespace std;

class Usr_Interface
{
        public:
                Usr_Interface () {}
                void cout_str(const char msg[]);
                void cout_str(string str);
                
                void cout_int(int x);
                void cout_float(float x);

                void new_line();
                
                void cout_str_no_space(string str);
                void cout_str_no_endl(const char msg[]);
                
                void cout_err(const char *err[]);
                void cout_err(const char err[]);
                string cin_str();
        private:
                int i;
};

#endif 
