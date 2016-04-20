/* This is console out stream handler implementation script 
 *
 * Linked to usr_interface.h
*/


// ======== includes ========

#include "usr_interface.h"

#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>

// ======= implementation ========

using namespace std;


// cout 

void Usr_Interface::new_line()
{
        cout<< endl;
}

void Usr_Interface::space()
{
        cout << endl << endl;
}

void Usr_Interface::cout_str(const char msg[])
{
        cout << msg << endl << endl;
}

void Usr_Interface::cout_str(string str)
{
        cout << str << endl << endl;
}

void Usr_Interface::cout_str_no_space(string str)
{
        cout << str << endl;
}

void Usr_Interface::cout_str_no_endl(const char msg[])
{
        cout << msg;
}

void Usr_Interface::cout_int(int x)
{
        cout << x << endl << endl;
}

void Usr_Interface::cout_float(float x)
{
        cout << x << endl << endl;
}

// error

void Usr_Interface::cout_err(const char *err[])
{
        cout << "An error occured in " << *(err)  << endl << endl << "------ Now exiting program ------" << endl << endl; 
}

void Usr_Interface::cout_err(const char err[])
{
        cout << "An error occured in " << err  << endl << endl << "------ Now exiting program ------" << endl << endl; 
}


//cin

string Usr_Interface::cin_str()
{
        string temp;
        cin >> temp;
        return temp;
}


//refresh

void Usr_Interface::start_DES_chrono()
{
         this->DES_start = chrono::high_resolution_clock::now();
}

int Usr_Interface::measure_DES_chrono()
{
         return floor( chrono::duration<double>( chrono::high_resolution_clock::now()-this->DES_start ).count() ) ; 
}

void Usr_Interface::refresh_DES_progress(int i, int max)
{
        int t = this->measure_DES_chrono();
        int p = ceil( (float) 100 * i / max);
        cout << flush << "\r" 
            << "DES->Progress: " << setw(3) << p << " \% -- " 
            << "Running for: " << setw(5) << t << "s -- "
            << "Predicted duration: " << setw(5) << ceil ( (float) 100 * t / p ) << "s";
}
