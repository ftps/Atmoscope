#include "tests.hpp"

int main(int argc, char* argv[])
{
    std::cout << "Test" << std::endl;
    std::cout << sizeof(double) << ", " << sizeof(long double) << std::endl;
    turbAtmos::testFlashMap();

    return 0;
}