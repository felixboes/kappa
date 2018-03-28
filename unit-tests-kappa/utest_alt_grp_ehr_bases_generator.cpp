#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "kappa/alt_grp_tuple.hpp"
#include "kappa/alt_grp_ehr_bases_generator.hpp"
#include <boost/math/special_functions/binomial.hpp>

BOOST_AUTO_TEST_SUITE(test_suit_alt_grp_ehr_bases_generator)


    BOOST_AUTO_TEST_CASE(test_constructor)
    {
        AltGrpTuple::parallel_case();
        AltGrpEhrBasesGenerator generator_par(5,1);                                 // g=5, m=1, h=11
        BOOST_CHECK(generator_par.h == 11);
        BOOST_CHECK(generator_par.possible_norm2permutations.size() == 23);         // = 2*h+1
        BOOST_CHECK(generator_par.bases.size() == 0);

        AltGrpTuple::radial_case();
        AltGrpEhrBasesGenerator generator_rad(5,1);                                 // g=5, m=1, h=10
        BOOST_CHECK(generator_rad.h == 10);
        BOOST_CHECK(generator_rad.possible_norm2permutations.size() == 21);         // = 2*h+1
        BOOST_CHECK(generator_rad.bases.size() == 0);
    }

    BOOST_AUTO_TEST_CASE(test_possible_norm2permutations)
    {
        AltGrpTuple::radial_case();
        AltGrpEhrBasesGenerator generator_rad(5,1);     //g=5, m=1, h=10
        generator_rad.generate_possible_norm2permutations();
        auto perms_rad = generator_rad.possible_norm2permutations;
        BOOST_CHECK(perms_rad.size() == 21);            // = 2*h +1
        BOOST_CHECK(perms_rad[0].size() == 0);          // expect perms_rad[ht].size() = 2*(ht choose 2) + 3*(ht choose 3)
        BOOST_CHECK(perms_rad[1].size() == 0);
        BOOST_CHECK(perms_rad[2].size() == 2);
        BOOST_CHECK(perms_rad[3].size() == 9);
        BOOST_CHECK(perms_rad[4].size() == 24);
        BOOST_CHECK(perms_rad[10].size() == 450);
        BOOST_CHECK(perms_rad[13].size() == 1014);
        BOOST_CHECK(perms_rad[16].size() == 1920);
        BOOST_CHECK(perms_rad[19].size() == 3249);
        //std::cout << "norm2permutations with height 4:" << std::endl;
        for(auto tau_rad : perms_rad[4])
        {
            BOOST_CHECK(PermutationManager::is_in_norm2notation(tau_rad));
            //std::cout << tau;
        }


        AltGrpTuple::parallel_case();
        AltGrpEhrBasesGenerator generator_par(5,1);     // g=5, m=1, h=11
        generator_par.generate_possible_norm2permutations();
        auto perms_par = generator_par.possible_norm2permutations;
        BOOST_CHECK(perms_par.size() == 23);            // = 2*h + 1
        BOOST_CHECK(perms_par[0].size() == 0);          // expect perms_par[ht].size() = 2*(ht-1 choose 2) + 3*(ht-1 choose 3)
        BOOST_CHECK(perms_par[1].size() == 0);
        BOOST_CHECK(perms_par[2].size() == 0);
        BOOST_CHECK(perms_par[3].size() == 2);
        BOOST_CHECK(perms_par[4].size() == 9);
        BOOST_CHECK(perms_par[10].size() == 324);
        BOOST_CHECK(perms_par[13].size() == 792);
        BOOST_CHECK(perms_par[16].size() == 1575);
        BOOST_CHECK(perms_par[19].size() == 2754);
        BOOST_CHECK(perms_par[22].size() == 4410);
        for(auto tau_par : perms_rad[4])
        {
            BOOST_CHECK(PermutationManager::is_in_norm2notation(tau_par));
        }
    }

    BOOST_AUTO_TEST_CASE(test_generate_bases_radial) {
        AltGrpTuple::radial_case();
        AltGrpEhrBasesGenerator generator_rad(1, 1);     // g=1, m=1, h=2
        auto bases_rad = generator_rad.generate_bases();
        BOOST_CHECK(bases_rad.size() == 3);
        BOOST_CHECK(bases_rad[0].size() == 0);
        BOOST_CHECK(bases_rad[1].size() == 0);
        BOOST_CHECK(bases_rad[2].size() == 1);
        BOOST_CHECK(bases_rad[3].size() == 2);
        BOOST_CHECK(bases_rad[4].size() == 1);
        BOOST_CHECK(bases_rad[2].basis.count( create_alt_grp_tuple(1, 1, 0, 2, 1) ) == 1);
        BOOST_CHECK(bases_rad[3].basis.count( create_alt_grp_tuple(1, 2, 1, 3, 2) ) == 1);
        BOOST_CHECK(bases_rad[3].basis.count( create_alt_grp_tuple(1, 2, 0, 3, 1) ) == 1);
        BOOST_CHECK(bases_rad[4].basis.count( create_alt_grp_tuple(1, 3, 1, 4, 2) ) == 1);
        //generator_rad.print_bases();


        AltGrpEhrBasesGenerator generator_rad2(1, 3);     // g=1, m=3, h=4
        auto bases_rad2 = generator_rad2.generate_bases();
        BOOST_CHECK(bases_rad2.size() == 7);
        BOOST_CHECK(bases_rad2[0].size() == 0);
        BOOST_CHECK(bases_rad2[1].size() == 0);
        BOOST_CHECK(bases_rad2[2].size() == 1);
        BOOST_CHECK(bases_rad2[3].size() == 57);
        BOOST_CHECK(bases_rad2[4].size() == 451);
        BOOST_CHECK(bases_rad2[5].size() == 1365);
        BOOST_CHECK(bases_rad2[6].size() == 1950);
        BOOST_CHECK(bases_rad2[7].size() == 1330);
        BOOST_CHECK(bases_rad2[8].size() == 350);
        BOOST_CHECK(bases_rad2[2].basis.count( create_alt_grp_tuple(2, 1, 0, 2, 1, 1, 0, 2, 1) ) == 1);
        for(uint8_t ht=0; ht<bases_rad2.size();ht++)
        {
            for(auto tuple : bases_rad2[ht].basis)
            {
                for(uint8_t i=1; i<=tuple.num_entries(); i++)
                {
                    BOOST_CHECK(PermutationManager::is_in_norm2notation(tuple.at(i)));
                }
            }
        }
        //generator_rad2.print_bases();
    }

    BOOST_AUTO_TEST_CASE(test_generate_bases_parallel)
    {
        AltGrpTuple::parallel_case();
        AltGrpEhrBasesGenerator generator_par(1,2);     // g=1, m=2, h=4
        auto bases_par = generator_par.generate_bases();
        BOOST_CHECK(bases_par.size() == 6);
        BOOST_CHECK(bases_par[0].size()==0);
        BOOST_CHECK(bases_par[1].size()==0);
        BOOST_CHECK(bases_par[2].size()==0);
        BOOST_CHECK(bases_par[3].size()==1);
        BOOST_CHECK(bases_par[4].size()==56);
        BOOST_CHECK(bases_par[5].size()==395);
        BOOST_CHECK(bases_par[6].size()==970);
        BOOST_CHECK(bases_par[7].size()==980);
        BOOST_CHECK(bases_par[8].size()==350);
        BOOST_CHECK(bases_par[3].basis.count( create_alt_grp_tuple(2, 2,1,3,2, 2,1,3,2) ) == 1);
        for(uint8_t ht=0; ht<bases_par.size();ht++)
        {
            for(auto tuple : bases_par[ht].basis)
            {
                for(uint8_t i=1; i<=tuple.num_entries(); i++)
                {
                    BOOST_CHECK(PermutationManager::is_in_norm2notation(tuple.at(i)));
                }
            }
        }
        //generator_par.print_bases();
    }



BOOST_AUTO_TEST_SUITE_END()