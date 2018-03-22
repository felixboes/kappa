#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "../kappa/permutation.hpp"

BOOST_AUTO_TEST_SUITE(test_suit_permutation)

    BOOST_AUTO_TEST_CASE(test_cycle_decomposition) {
        Permutation tau(6);
        tau[0] = 0;
        tau[1] = 2;
        tau[2] = 4;
        tau[3] = 5;
        tau[4] = 1;
        tau[5] = 3;

        auto cycle_decomp = tau.cycle_decomposition();

        BOOST_CHECK(tau.size() == 6);
        BOOST_CHECK(cycle_decomp.size() == 3);
        BOOST_CHECK(cycle_decomp.count(0) == 1);
        BOOST_CHECK(cycle_decomp.count(1) == 1);
        BOOST_CHECK(cycle_decomp.count(3) == 1);

        Permutation cycle0 = cycle_decomp.at(0);
        BOOST_CHECK(cycle0.at(0) == 0 and cycle0.at(1) == 6 and cycle0.at(2) == 6 and cycle0.at(3) == 6
                    and cycle0.at(4) == 6 and cycle0.at(5) == 6);

        Permutation cycle1 = cycle_decomp.at(1);
        BOOST_CHECK(cycle1.at(0) == 6 and cycle1.at(1) == 2 and cycle1.at(2) == 4 and cycle1.at(3) == 6
                    and cycle1.at(4) == 1 and cycle1.at(5) == 6);

        Permutation cycle3 = cycle_decomp.at(3);
        BOOST_CHECK(cycle3.at(0) == 6 and cycle3.at(1) == 6 and cycle3.at(2) == 6 and cycle3.at(3) == 5
                    and cycle3.at(4) == 6 and cycle3.at(5) == 3);
    }

    BOOST_AUTO_TEST_CASE(test_long_cycle) {
        Permutation long_cycle = Permutation::long_cycle(15);
        Permutation long_cycle_inv = Permutation::long_cycle_inv(15);

        BOOST_CHECK(long_cycle.size() == 16);

        for (int i = 0; i < 15; i++) {
            BOOST_CHECK(long_cycle[i] == i + 1);
            BOOST_CHECK(long_cycle_inv[i + 1] == i);
        }
        BOOST_CHECK(long_cycle[15] == 0);
        BOOST_CHECK(long_cycle_inv[0] == 15);
    }

    BOOST_AUTO_TEST_SUITE_END()