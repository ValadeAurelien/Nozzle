/* This is console out stream handler implementation script 
 *
 * Linked to usr_interface.h
*/


// ======== includes ========

#include "../include/usr_interface.h"

#include <iostream>
#include <iomanip>

// ======= implementation ========

using namespace std;

void Usr_Interface::cout_str(const char msg[])
{
        cout << msg << endl << endl;
}

void Usr_Interface::cout_str(string str)
{
        cout << str << endl << endl;
}

void Usr_Interface::new_line()
{
        cout<< endl;
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

void Usr_Interface::cout_err(const char *err[])
{
        cout << "An error occured in " << *(err)  << endl << endl << "------ Now exiting program ------" << endl << endl; 
}

void Usr_Interface::cout_err(const char err[])
{
        cout << "An error occured in " << err  << endl << endl << "------ Now exiting program ------" << endl << endl; 
}

string Usr_Interface::cin_str()
{
        string temp;
        cin >> temp;
        return temp;
}
