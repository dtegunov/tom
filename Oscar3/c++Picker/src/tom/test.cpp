



#include <helper/filesystem.hpp>


#include <iostream>


int main(int argc, char **argv) {



    helper::fs::path p("/afs/ipp/home/t/thaller/usr/src/boost_1_34_1");

    std::cout << p << ": " << helper::fs::is_directory(p) << std::endl;
    std::cout << (p/"rst.css") << ": " << helper::fs::is_directory((p/"rst.css")) << " " << helper::fs::is_regular((p/"rst.css")) << std::endl;



    std::cout << helper::fs::create_directory("TEST1") << std::endl;


    return 0;

}





