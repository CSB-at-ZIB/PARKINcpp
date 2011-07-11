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
        case '1' : testlinalg();        break;
        case '2' : testnonlin();        break;
        case '3' : testsystem();        break;
        case '4' : testsystem_aux();    break;
        case '5' : testpfizer_simple(); break;
        case '6' :
        default  : testparkin_aux();    break;
    }

    return 0;
}
