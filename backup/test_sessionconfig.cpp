#include <iostream>

#include <boost/program_options.hpp>

void print_usage(int argc, char** argv, boost::program_options::options_description& desc)
{
    std::cout << "Usage: " << argv[0] << " -g arg -m arg optional_args" << std::endl;
    std::cout << desc << std::endl;
}

int main( int argc, char** argv )
{
    // Declare the supboost::program_optionsrted options.
    boost::program_options::options_description desc("Command line options");
    desc.add_options()
            ("genus,g", boost::program_options::value<int>(), "the genus of the Riemannian surfaces")
            ("punctures,m", boost::program_options::value<int>(), "the number of punctures of the Riemannian surfaces")
            ("rat,r", "uses rational numbers as coefficients")
            ("fin,n", boost::program_options::value<int>(), "uses the finite field F_n with n a prime number")
            ("num_threads,t", boost::program_options::value<int>(), "the number of threads used in matrix computations")
            ("first_diff", boost::program_options::value<int>(), "start with the differential first_diff")
            ("last_diff", boost::program_options::value<int>(), "end with the differential last_diff")
            ("help,h", "produce help message")
    ;

    boost::program_options::variables_map vm;

    try{
        boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);

        if (vm.count("help")) {
            print_usage(argc, argv, desc);
            return 0;
        }

        try{
            boost::program_options::notify(vm);
        }
        catch(boost::program_options::required_option& e)
        {
            std::cerr << "Error: " << e.what() << std::endl << std::endl;
            print_usage(argc, argv, desc);
            return 1;
        }
    }
    catch(boost::program_options::error& e)
    {
        std::cerr << "Error: " << e.what() << std::endl << std::endl;
        print_usage(argc, argv, desc);
        return 1;
    }

    std::cout << "genus = " << vm["genus"].as<int>() << "; punctures = " << vm["punctures"].as<int>() << std::endl;

    return 0;
}
