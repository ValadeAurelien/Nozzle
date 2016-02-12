/* This is console out stream handler implementation script 
 *
 * Linked to usr_interface.h
*/


// ======== includes ========

#include "usr_interface.h"


// ======= implementation ========

using namespace std;

Usr_Interface::cout_str(string msg)
{
        cout << msg << endl;
}

string Usr_Interface::cin_str()
{
        string temp;
        cin >> temp;
        return temp;
}
