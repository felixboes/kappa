#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "../kappa/permutation_manager.hpp"

BOOST_AUTO_TEST_SUITE(test_suit_permutation_manager)

    BOOST_AUTO_TEST_CASE(test_disjoint_transpositions) {
        Transposition t1(1, 1), t2(1, 2), t3(2, 5), t4(3, 4);

        BOOST_CHECK(!PermutationManager::are_symbolwise_disjoint(t1, t2));
        BOOST_CHECK(PermutationManager::are_symbolwise_disjoint(t1, t4));
        BOOST_CHECK(PermutationManager::are_symbolwise_disjoint(t2, t4));
        BOOST_CHECK(!PermutationManager::are_symbolwise_disjoint(t3, t2));


        BOOST_CHECK(PermutationManager::have_disjoint_support(t1, t2));
        BOOST_CHECK(PermutationManager::have_disjoint_support(t1, t4));
        BOOST_CHECK(PermutationManager::have_disjoint_support(t2, t4));
        BOOST_CHECK(!PermutationManager::have_disjoint_support(t3, t2));
    }

    BOOST_AUTO_TEST_CASE(test_transposition_value_at) {
        Transposition t1(1, 1), t2(1, 2);

        BOOST_CHECK(PermutationManager::value_at(t1, 1) == 1);
        BOOST_CHECK(PermutationManager::value_at(t1, 5) == 5);

        BOOST_CHECK(PermutationManager::value_at(t2, 1) == 2);
        BOOST_CHECK(PermutationManager::value_at(t2, 2) == 1);
        BOOST_CHECK(PermutationManager::value_at(t2, 3) == 3);
    }

    BOOST_AUTO_TEST_CASE(test_norm2permutation_value_at) {
        Norm2Permutation t1(Transposition(1, 1), Transposition(2, 3));
        Norm2Permutation t2(Transposition(1, 2), Transposition(2, 3));

        BOOST_CHECK(PermutationManager::value_at(t1, 1) == 1);
        BOOST_CHECK(PermutationManager::value_at(t1, 2) == 3);
        BOOST_CHECK(PermutationManager::value_at(t1, 3) == 2);
        BOOST_CHECK(PermutationManager::value_at(t1, 5) == 5);

        BOOST_CHECK(PermutationManager::value_at(t2, 1) == 2);
        BOOST_CHECK(PermutationManager::value_at(t2, 2) == 3);
        BOOST_CHECK(PermutationManager::value_at(t2, 3) == 1);
        BOOST_CHECK(PermutationManager::value_at(t2, 4) == 4);
    }

    BOOST_AUTO_TEST_CASE(test_num_cycles) {
        Permutation tau(6);
        tau[0] = 0;
        tau[1] = 2;
        tau[2] = 4;
        tau[3] = 5;
        tau[4] = 1;
        tau[5] = 3;
        BOOST_CHECK(PermutationManager::num_cycles(tau) == 3);
    }

    BOOST_AUTO_TEST_CASE(test_descending_cycle_decomposition_of_product) {
        Norm2Permutation t1(Transposition(1, 4), Transposition(4, 0));
        Norm2Permutation t2(Transposition(0, 4), Transposition(3, 5));
        Norm2Permutation t3(Transposition(4, 3), Transposition(1, 6));

        auto cycle_decomp_t1_t2 = PermutationManager::descending_cycle_decomposition_of_product(t1, t2);
        BOOST_CHECK(cycle_decomp_t1_t2.size() == 2);

        auto cycle0 = cycle_decomp_t1_t2[0];
        BOOST_CHECK(cycle0.size() == 2 and cycle0[0] == 5 and cycle0[1] == 3);

        auto cycle1 = cycle_decomp_t1_t2[1];
        BOOST_CHECK(cycle1.size() == 2 and cycle1[0] == 4 and cycle1[1] == 1);


        auto cycle_decomp_t2_t3 = PermutationManager::descending_cycle_decomposition_of_product(t2, t3);
        BOOST_CHECK(cycle_decomp_t2_t3.size() == 2);

        auto cycle2 = cycle_decomp_t2_t3[0];
        BOOST_CHECK(cycle2.size() == 2 and cycle2[0] == 6 and cycle2[1] == 1);

        auto cycle3 = cycle_decomp_t2_t3[1];
        BOOST_CHECK(cycle3.size() == 4 and cycle3[0] == 5 and cycle3[1] == 3 and cycle3[2] == 0 and cycle3[3] == 4);
    }

    BOOST_AUTO_TEST_CASE(test_inverse_permutation) {
        Permutation tau(6);
        tau[0] = 0;
        tau[1] = 2;
        tau[2] = 4;
        tau[3] = 5;
        tau[4] = 1;
        tau[5] = 3;
        Permutation inv = PermutationManager::inverse(tau);

        BOOST_CHECK(inv.size() == 6);
        BOOST_CHECK(inv[0] == 0 and inv[1] == 4 and inv[2] == 1 and inv[3] == 5 and inv[4] == 2 and inv[5] == 3);
    }

    BOOST_AUTO_TEST_CASE(test_multiply_with_permutation_manager) {
        Permutation tau(6);
        tau[0] = 0;
        tau[1] = 2;
        tau[2] = 4;
        tau[3] = 5;
        tau[4] = 1;
        tau[5] = 3;

        Transposition t1(2, 5), t2(2, 0), t3(3, 3);
        Norm2Permutation sigma(t1, t2);

        Permutation tau_t1(tau);
        PermutationManager::multiply(tau_t1, t1);
        BOOST_CHECK(tau_t1.size() == 6);
        BOOST_CHECK(tau_t1[0] == 0 and tau_t1[1] == 2 and tau_t1[2] == 3 and tau_t1[3] == 5 and tau_t1[4] == 1
                    and tau_t1[5] == 4);

        Permutation tau_t2(tau);
        PermutationManager::multiply(tau_t2, t2);
        BOOST_CHECK(tau_t2.size() == 6);
        BOOST_CHECK(tau_t2[0] == 4 and tau_t2[1] == 2 and tau_t2[2] == 0 and tau_t2[3] == 5 and tau_t2[4] == 1
                    and tau_t2[5] == 3);

        Permutation tau_t3(tau);
        PermutationManager::multiply(tau_t3, t3);
        BOOST_CHECK(tau_t3 == tau);

        Permutation tau_sigma(tau);
        PermutationManager::multiply(tau_sigma, sigma);
        BOOST_CHECK(tau_t2.size() == 6);
        BOOST_CHECK(tau_sigma[0] == 3 and tau_sigma[1] == 2 and tau_sigma[2] == 0 and tau_sigma[3] == 5
                    and tau_sigma[4] == 1 and tau_sigma[5] == 4);
    }

    BOOST_AUTO_TEST_CASE(test_norm2notation_of_norm2permutation)
    {
        Norm2Permutation t(Transposition(3,1), Transposition(1,2));
        BOOST_CHECK(PermutationManager::has_norm_two(t));
        BOOST_CHECK(!PermutationManager::is_in_norm2notation(t));
        BOOST_CHECK(PermutationManager::write_in_norm2notation(t));
        BOOST_CHECK(PermutationManager::is_in_norm2notation(t));
        BOOST_CHECK(t == Norm2Permutation(Transposition(2,1), Transposition(3,2)));

        Norm2Permutation s1(Transposition(3,1), Transposition(2,2));
        BOOST_CHECK(!PermutationManager::has_norm_two(s1));
        BOOST_CHECK(!PermutationManager::is_in_norm2notation(s1));
        BOOST_CHECK(!PermutationManager::write_in_norm2notation(s1));

        Norm2Permutation s2(Transposition(3,1), Transposition(1,3));
        BOOST_CHECK(!PermutationManager::has_norm_two(s2));
        BOOST_CHECK(!PermutationManager::is_in_norm2notation(s2));
        BOOST_CHECK(!PermutationManager::write_in_norm2notation(s2));

        Norm2Permutation s3(Transposition(3,1), Transposition(3,1));
        BOOST_CHECK(!PermutationManager::has_norm_two(s3));
        BOOST_CHECK(!PermutationManager::is_in_norm2notation(s3));
        BOOST_CHECK(!PermutationManager::write_in_norm2notation(s3));
    }

    BOOST_AUTO_TEST_CASE(test_substitute_index_of_permutations)
    {
        Transposition t(2,17);
        PermutationManager::substitute_index(t, 2, 3);
        BOOST_CHECK(t == Transposition(3,17));
        PermutationManager::substitute_index(t, 5, 3);
        BOOST_CHECK(t == Transposition(3,17));

        Norm2Permutation s(Transposition(3,1), Transposition(1,2));
        PermutationManager::substitute_index(s, 3,15);
        BOOST_CHECK(s == Norm2Permutation(Transposition(15,1), Transposition(1,2)));
        PermutationManager::substitute_index(s, 1,2);
        BOOST_CHECK(s == Norm2Permutation(Transposition(15,2), Transposition(2,2)));
        PermutationManager::substitute_index(s, 15,15);
        BOOST_CHECK(s == Norm2Permutation(Transposition(15,2), Transposition(2,2)));
    }

    BOOST_AUTO_TEST_CASE(test_drop_fixed_point)
    {
        Transposition t(2,17);
        PermutationManager::drop_fixed_point(t, 3);
        BOOST_CHECK(t == Transposition(2,16));
        PermutationManager::drop_fixed_point(t, 18);
        BOOST_CHECK(t == Transposition(2,16));
        PermutationManager::drop_fixed_point(t, 0);
        BOOST_CHECK(t == Transposition(1,15));

        Norm2Permutation s(Transposition(3,8), Transposition(1,10));
        PermutationManager::drop_fixed_point(s, 15);
        BOOST_CHECK(s == Norm2Permutation(Transposition(3,8), Transposition(1,10)));
        PermutationManager::drop_fixed_point(s, 7);
        BOOST_CHECK(s == Norm2Permutation(Transposition(3,7), Transposition(1,9)));
    }

    BOOST_AUTO_TEST_CASE(test_orientation_sign)
    {
        Permutation sigma(5);
        sigma[0] = 4;
        sigma[4] = 0;
        sigma[1] = 3;
        sigma[3] = 1;
        sigma[2] = 2;
        std::map< uint8_t, int8_t > orientation_sign = PermutationManager::orientation_sign_for_ehrenfried(sigma);

        BOOST_CHECK(orientation_sign.size()==5);
        BOOST_CHECK(orientation_sign.at(0)==1);
        BOOST_CHECK(orientation_sign.at(1)==-1);
        BOOST_CHECK(orientation_sign.at(2)==0);
        BOOST_CHECK(orientation_sign.at(3)==1);
        BOOST_CHECK(orientation_sign.at(4)==1);
    }



BOOST_AUTO_TEST_SUITE_END()