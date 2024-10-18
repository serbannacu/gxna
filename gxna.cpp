#include "Args.h"
#include "Experiment.h"
#include "Exception.h"

#include <fstream>
#include <iostream>

int main(int argc, char *argv[]) {
    gxna::Args args;
    try {
        args.read("gxna.arg", false); // no error if file does not exist
        args.parse(argc, argv);
        args.check();
    }
    catch (const std::exception& e) {
        std::cerr << "Args error: " << e.what() << std::endl;
        return 1;
    }

    std::cerr << "GXNA Version 3.0\n";
    srand(args.seed);

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
