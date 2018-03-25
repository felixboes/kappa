#include "alt_grp_ehr_bases_generator.hpp"

AltGrpEhrBasesGenerator::AltGrpEhrBasesGenerator(uint8_t _g, uint8_t _m) :
    g(_g),
    m(_m),
    h(2*_g+_m-1+AltGrpTuple::get_min_symbol()),
    possible_norm2permutations(2*h+1, std::vector< Norm2Permutation >()),
    bases(2*h+1, std::vector< AltGrpTuple >())
{
    h=2*g+m;
    if(AltGrpTuple::is_radial())
    {
        h--;
    }
}

std::vector< std::vector< AltGrpTuple >>& AltGrpEhrBasesGenerator::generate_bases()
{
    if(h % 2)
    {
        std::cerr << "Error in std::vector< std::vector< AltGrpTuple >>& AltGrpEhrBasesGenerator::generate_bases(): h="
                 << std::to_string(h) << " is odd." << std::endl;
        return bases;
    }
    if( (possible_norm2permutations.size() != (unsigned)2*h+1) or (bases.size() != (unsigned)2*h+1) )
    {
        std::cerr << "Error in std::vector< std::vector< AltGrpTuple >>& AltGrpEhrBasesGenerator::generate_bases(): "
                 << " possible_norm2permutations.size()=" << std::to_string(possible_norm2permutations.size())
                 << " and bases.size()= " << std::to_string(bases.size()) << " but both should equal 2*h + 1 = "
                 << std::to_string(2*h+1) << std::endl;
        return bases;
    }
    AltGrpTuple empty_tuple(static_cast<size_t>(h / 2));
    generate_possible_norm2permutations();
    generate_bases_recursively(0, 0, empty_tuple);
    return bases;
}

void AltGrpEhrBasesGenerator::generate_possible_norm2permutations()
{
    uint8_t min_symbol = AltGrpTuple::get_min_symbol();
    Norm2Permutation tau;

    //compute all 3-cycles
    for(uint8_t height = 2+min_symbol; height <= 2*h; height++)
    {
        for(uint8_t a = 1+min_symbol; a < height; a++)
        {
            for(uint8_t b = min_symbol; b < a; b++)
            {
                //two possible 3-cycles with symbols height > a > b
                tau.first = Transposition(a,b);
                tau.second = Transposition(height,a);
                possible_norm2permutations[height].push_back(tau);  // (height b a) = (a b)(height a)

                tau.second = Transposition(height,b);
                possible_norm2permutations[height].push_back(tau);  // (height a b) = (a b)(height b)
            }
        }
    }

    //compute all products of disjoint transpositions
    for(uint8_t height = 3+min_symbol; height <= 2*h; height++)
    {
        for(uint8_t a = 2+min_symbol; a < height; a++)
        {
            for(uint8_t b = 1+min_symbol; b < a; b++)
            {
                for(uint8_t c = min_symbol; c < b; c++)
                {
                    //three possible norm2permutations with symbols height > a > b > c
                    tau.first = Transposition(a,b);
                    tau.second = Transposition(height,c);
                    possible_norm2permutations[height].push_back(tau);  // (a b)(height c)

                    tau.first = Transposition(a,c);
                    tau.second = Transposition(height,b);
                    possible_norm2permutations[height].push_back(tau);  // (a c)(height b)

                    tau.first = Transposition(b,c);
                    tau.second = Transposition(height,a);
                    possible_norm2permutations[height].push_back(tau);  // (b c)(height a)
                }
            }
        }
    }
}

void AltGrpEhrBasesGenerator::generate_bases_recursively(uint8_t curr_num_entries, uint8_t curr_max_symbol,
                                                        AltGrpTuple& curr_tuple)
{
    if(curr_num_entries < h/2)
    {
        uint8_t new_num_entries = curr_num_entries +1;
        uint8_t new_max_symbol;

        // the height of the new entry must be >= the first symbol of the current entry in norm2notation, as we require
        // unstablity
        uint8_t new_min_height = 0;
        if(curr_num_entries > 0)
        {
            new_min_height = curr_tuple.at(curr_num_entries).first.first;
        }

        // iterate with all possibilities for the new entry
        for(uint8_t height= new_min_height; height<=2*h; height++)
        {
            for (const Norm2Permutation &tau : possible_norm2permutations[height])
            {
                new_max_symbol = std::max(curr_max_symbol, height);
                curr_tuple.p =new_max_symbol;
                curr_tuple.at(new_num_entries) = tau;
                generate_bases_recursively(new_num_entries, new_max_symbol, curr_tuple);
            }
        }
    }
    else if(curr_num_entries ==h/2)
    {
        //if curr_tuple is a basis element, add it to the corresponding basis
        if( curr_tuple.has_correct_num_cycles(m)
            and curr_tuple.has_no_common_fixed_points_except_zero())
        {
            bases.at(curr_max_symbol).push_back(curr_tuple);
        }
    }
}

void AltGrpEhrBasesGenerator::print_bases()
{
    std::cout << "Bases of th Ehrenfried complex for g=" << std::to_string(g) << ", m=" << std::to_string(m) << ", h="
              << std::to_string(h) << "in the " << ( AltGrpTuple::is_radial() ? "radial" : "parallel") << " case."
              << std::endl;
    for(uint8_t p=0; p < bases.size(); p++)
    {
        std::cout << std::endl << "The " << std::to_string(p) << "-th module has dimension "
                  << std::to_string(bases[p].size()) << " and the following basis elements:" << std::endl;
        for(auto tuple : bases[p])
        {
            std::cout << tuple << std::endl;
        }
    }
}

