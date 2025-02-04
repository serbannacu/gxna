#include "Args.h"
#include "Experiment.h"
#include "Exception.h"

#include <fstream>
#include <iostream>

int main(int argc, char *argv[]) {
    std::cout << "GXNA Version 3.0" << std::endl;

    gxna::Args args;
    try {
        args.parse(argc, argv);
        args.check();
    }
    catch (const std::exception& e) {
        std::cerr << "Args error: " << e.what() << std::endl;
        args.usage(std::cerr);
        return 1;
    }

    try {
        gxna::Experiment experiment(args);
        experiment.run();
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    return 0;
}
