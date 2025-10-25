#include "../core/matrix.h"
#include <iostream>
using namespace  gons::utilites::LOG_MSG;



int main(int argc, char **argv)
{
    UNUSED(argc);
    UNUSED(argv);

    gons::Matrix<double, 2, 2> m1 = {
        {1, 2},
        {3, 4}
    };
    
    gons::Matrix<double, 2, 2> m2 = m1.Inverse();
    m2.Print("m2 Inverse : \n");

    gons::Matrix<double, 3, 3> m3 = {
        {1, 2, 3},  
        {4, 5, 6},
        {7, 8, 9}
    };
    
    gons::Matrix<double, 3, 3> m4 = m3.Inverse();
    m4.Print("m4 Inverse : \n");
    
    return 0;

}
