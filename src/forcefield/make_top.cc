#include "ff_lib.cpp"
#include <iostream>
#include <fstream>
#include <string>

void print_usage ()
{
    std::string line;
    std::string fn;
    fn = std::string(__CURRDIR__) + "/" + std::string(__FILE__) + ".usage";
    std::ifstream usg_fp (fn);
    while (std::getline(usg_fp, line))
        std::cout << line << std::endl;
    usg_fp.close();
    return;
}

int main (int argc, char** argv)
{
    print_usage();
    return 0;
}
