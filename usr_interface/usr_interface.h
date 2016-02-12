
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

#include <iostream>
#include <iomanip>
#include <string>

// ======== interface ========

class Usr_Interface
{
        public:
                void cout_str(string msg);
                string cin_str();
        private:
                int i;
};

#endif 
