#include <iostream>
#include "src/random.h"

int main(int argc, char **argv)
{
    Crandom ran(1);
    for (int i = 0; i < 10; i++)
    {
        std::cout << ran.random() << std::endl;
    }
    return 0;
}