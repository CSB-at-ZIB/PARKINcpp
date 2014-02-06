/// a comment for documentation

#include <iostream>
#include "testcases.h"

using namespace std;

int main(int argc, char* argv[])
{
    cout << "Hello world!" << endl;

    ///

    char ch = (argc > 1) ? *argv[1] : ' ';

    switch(ch)
    {
        case '0' : testmshoot();            break;
        case '1' : testlinalg();            break;
        case '2' : testnonlin();            break;
        case '3' : testsystem();            break;
        case '4' : testsystem_aux();        break;
        case '5' : testpfizer_simple();     break;
        //case '6' : testfoerster_react_b();  break;
        case '6' : testfoerster_react_c();  break;
        case '7' :
        default  : testparkin_aux();    break;
    }

    return 0;
}
