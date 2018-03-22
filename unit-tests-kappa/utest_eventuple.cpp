#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "../kappa/eventuple.hpp"

BOOST_AUTO_TEST_SUITE(test_suit_eventuple)

    class TestNorm2Permutations{
        public:
            TestNorm2Permutations() {}

            //provide several Norm2Permutations for the tests. Always denote by s_i the Norm2Notation of t_i.

            Norm2Permutation t0 = Norm2Permutation(Transposition(0,3),Transposition(2,3));
            Norm2Permutation s0 = Norm2Permutation(Transposition(2,0),Transposition(3,0));

            Norm2Permutation t1 = Norm2Permutation(Transposition(1,3),Transposition(2,3));
            Norm2Permutation s1 = Norm2Permutation(Transposition(2,1),Transposition(3,1));

            Norm2Permutation t2 = Norm2Permutation(Transposition(1,2),Transposition(3,4));
            Norm2Permutation s2 = Norm2Permutation(Transposition(2,1),Transposition(4,3));

            Norm2Permutation t3 = Norm2Permutation(Transposition(3,5),Transposition(6,1));
            Norm2Permutation s3 = Norm2Permutation(Transposition(5,3),Transposition(6,1));

            //t4 has norm zero and hence cannot be written in Norm2Notation
            Norm2Permutation t4 = Norm2Permutation(Transposition(3,5),Transposition(5,3));

            Norm2Permutation t5 = Norm2Permutation(Transposition(2,4),Transposition(4,5));
            Norm2Permutation s5 = Norm2Permutation(Transposition(4,2),Transposition(5,4));

            Norm2Permutation t6 = Norm2Permutation(Transposition(3,5),Transposition(6,0));
            Norm2Permutation s6 = Norm2Permutation(Transposition(5,3),Transposition(6,0));

            Norm2Permutation t7 = Norm2Permutation(Transposition(4,2),Transposition(7,8));
            Norm2Permutation s7 = Norm2Permutation(Transposition(4,2),Transposition(8,7));

            Norm2Permutation t8 = Norm2Permutation(Transposition(5,3),Transposition(7,8));
            Norm2Permutation s8 = Norm2Permutation(Transposition(5,3),Transposition(8,7));

            Norm2Permutation s9 = Norm2Permutation(Transposition(3,2),Transposition(4,2));
    };

    BOOST_FIXTURE_TEST_CASE(test_access_eventuples_and_check_norm, TestNorm2Permutations)
    {
        EvenTuple::parallel_case();
        BOOST_CHECK(!EvenTuple::is_radial());

        EvenTuple T(6,3);
        BOOST_CHECK(T.num_norm2permutations()==3);

        T.at(1)=t1;
        T.at(2)=t2;
        T[3]=t3;
        BOOST_CHECK(T.num_norm2permutations()==3);
        BOOST_CHECK(T.has_correct_norm());

        T.at(3)=t4;
        BOOST_CHECK(!T.has_correct_norm());
    }

    BOOST_FIXTURE_TEST_CASE(test_eventuple_norm2notation, TestNorm2Permutations)
    {
        EvenTuple T(6,3);
        T.at(1)=t1;
        T.at(2)=t2;
        T.at(3)=t3;

        BOOST_CHECK(T.write_in_Norm2Notation());
        BOOST_CHECK(T.at(1)== s1);
        BOOST_CHECK(T.at(2)== s2);
        BOOST_CHECK(T.at(3)== s3);

        T.at(3)=t4;
        BOOST_CHECK(!T.write_in_Norm2Notation());
    }

    BOOST_FIXTURE_TEST_CASE(test_monotone_eventuple_parallelcase, TestNorm2Permutations)
    {
        EvenTuple::parallel_case();
        EvenTuple S1(6,3), S2(6,3), S3(6,3);

        S1.at(1)=t1;
        S1.at(2)=t5;
        S1.at(3)=t3;
        BOOST_CHECK(S1.monotone());

        S2.at(1)=t1;
        S2.at(2)=t3;
        S2.at(3)=t5;
        BOOST_CHECK(S2.monotone());

        S3.at(1)=t5;
        S3.at(2)=t1;
        S3.at(3)=t3;
        BOOST_CHECK(!S3.monotone());
    }

    BOOST_FIXTURE_TEST_CASE(test_monotone_eventuple_radialcase, TestNorm2Permutations)
    {
        EvenTuple::radial_case();
        EvenTuple S1(6,3), S2(6,3), S3(6,3);

        S1.at(1)=t0;
        S1.at(2)=t5;
        S1.at(3)=t6;
        BOOST_CHECK(S1.monotone());

        S2.at(1)=t0;
        S2.at(2)=t6;
        S2.at(3)=t5;
        BOOST_CHECK(S2.monotone());

        S3.at(1)=t5;
        S3.at(2)=t0;
        S3.at(3)=t6;
        BOOST_CHECK(not S3.monotone());
    }

    BOOST_FIXTURE_TEST_CASE(test_f_for_eventuples_parallelcase, TestNorm2Permutations)
    {
        EvenTuple::parallel_case();
        EvenTuple S1(6,3), S2(6,3), S3(8,2), S4(6,2);

        S1.at(1)=t1;
        S1.at(2)=t5;
        S1.at(3)=t3;
        BOOST_CHECK(S1.f(1));
        BOOST_CHECK( S1.at(1) == Norm2Permutation(Transposition(4,3),Transposition(5,4)) );
        BOOST_CHECK( S1.at(2) == t1 );
        BOOST_CHECK( S1.at(3) == t3 );

        S2.at(1)=t1;
        S2.at(2)=t5;
        S2.at(3)=t3;
        BOOST_CHECK(S2.f(2));
        BOOST_CHECK( S2.at(1) == t1 );
        BOOST_CHECK( S2.at(2) == s3);
        BOOST_CHECK( S2.at(3) == Norm2Permutation(Transposition(4,3),Transposition(2,3)) );

        S3.at(1)=t3;
        S3.at(2)=t7;
        BOOST_CHECK(S3.f(1));
        BOOST_CHECK( S3.at(1) == Norm2Permutation(Transposition(6,1),Transposition(8,7)) );
        BOOST_CHECK( S3.at(2) == Norm2Permutation(Transposition(4,2),Transposition(5,3)) );

        S4.at(1)=t7;
        S4.at(2)=t5;
        BOOST_CHECK(!S4.f(1));
    }

    BOOST_FIXTURE_TEST_CASE(test_f_for_eventuples_radialcase, TestNorm2Permutations)
    {
        EvenTuple::radial_case();
        BOOST_CHECK(EvenTuple::is_radial());
        EvenTuple S1(6,3), S2(6,3), S3(8,2), S4(6,2);
        Norm2Permutation t7(Transposition(4,2),Transposition(7,8));

        S1.at(1)=t0;
        S1.at(2)=t5;
        S1.at(3)=t6;
        BOOST_CHECK(S1.f(1));
        BOOST_CHECK( S1.at(1) == Norm2Permutation(Transposition(4,3),Transposition(5,4)) );
        BOOST_CHECK( S1.at(2) == t0);
        BOOST_CHECK( S1.at(3) == t6 );

        S2.at(1)=t0;
        S2.at(2)=t5;
        S2.at(3)=t6;
        BOOST_CHECK(S2.f(2));
        BOOST_CHECK( S2.at(1) == t0 );
        BOOST_CHECK( S2.at(2) == s6);
        BOOST_CHECK( S2.at(3) == Norm2Permutation(Transposition(4,3),Transposition(2,3)) );

        S3.at(1)=t6;
        S3.at(2)=t7;
        BOOST_CHECK(S3.f(1));
        BOOST_CHECK( S3.at(1) == Norm2Permutation(Transposition(6,0),Transposition(8,7)) );
        BOOST_CHECK( S3.at(2) == Norm2Permutation(Transposition(4,2),Transposition(5,3)) );

        S4.at(1)=t7;
        S4.at(2)=t5;
        BOOST_CHECK(!S4.f(1));
    }

    BOOST_FIXTURE_TEST_CASE(test_phi_for_eventuples_parallelcase, TestNorm2Permutations)
    {
        EvenTuple::parallel_case();
        EvenTuple S(6,3), S11(6,3), S21(6,3), S32(6,3), S13(6,3), T(6,3), N(6,3);

        S.at(1)=t1;
        S.at(2)=t5;
        S.at(3)=t3;
        N.at(1)=s1;
        N.at(2)=s5;
        N.at(3)=s3;

        S11=S;
        BOOST_CHECK(S11.phi(1,1));
        BOOST_CHECK(S11 == N);
        BOOST_CHECK(S11.phi(2,2));
        BOOST_CHECK(S11 == N);
        BOOST_CHECK(S11.phi(3,3));
        BOOST_CHECK(S11 == N);

        S21=S;
        BOOST_CHECK(S21.phi(2,1));
        BOOST_CHECK( S21.at(1) == Norm2Permutation(Transposition(4,3),Transposition(5,4)) );
        BOOST_CHECK( S21.at(2) == s1 );
        BOOST_CHECK( S21.at(3) == s3 );

        S32=S;
        BOOST_CHECK(S32.phi(3,2));
        BOOST_CHECK( S32.at(1) == s1 );
        BOOST_CHECK( S32.at(2) == s3 );
        BOOST_CHECK( S32.at(3) == s9 );

        S13=S;
        BOOST_CHECK(S13.phi(3,1));;
        BOOST_CHECK( S13.at(1) == Norm2Permutation(Transposition(5,1),Transposition(6,2)) );
        BOOST_CHECK( S13.at(2) == s1 );
        BOOST_CHECK( S13.at(3) == s9 );

        T.at(1)=t8;
        T.at(2)=t5;
        T.at(3)=t3;
        BOOST_CHECK(!T.phi(3,1));
    }

    BOOST_FIXTURE_TEST_CASE(test_phi_for_eventuples_radialcase, TestNorm2Permutations)
    {
        EvenTuple::radial_case();
        EvenTuple S(6,3), S11(6,3), S21(6,3), S32(6,3), S13(6,3), T(6,3), N(6,3);

        S.at(1)=t0;
        S.at(2)=t5;
        S.at(3)=t6;
        N.at(1)=s0;
        N.at(2)=s5;
        N.at(3)=s6;

        S11=S;
        BOOST_CHECK(S11.phi(1,1));
        BOOST_CHECK(S11 == N);
        BOOST_CHECK(S11.phi(2,2));
        BOOST_CHECK(S11 == N);
        BOOST_CHECK(S11.phi(3,3));
        BOOST_CHECK(S11 == N);

        S21=S;
        BOOST_CHECK(S21.phi(2,1));
        BOOST_CHECK( S21.at(1) == Norm2Permutation(Transposition(4,3),Transposition(5,4)) );
        BOOST_CHECK( S21.at(2) == s0 );
        BOOST_CHECK( S21.at(3) == s6 );

        S32=S;
        BOOST_CHECK(S32.phi(3,2));
        BOOST_CHECK( S32.at(1) == s0 );
        BOOST_CHECK( S32.at(2) == s6 );
        BOOST_CHECK( S32.at(3) == s9 );

        S13=S;
        BOOST_CHECK(S13.phi(3,1));;
        BOOST_CHECK( S13.at(1) == Norm2Permutation(Transposition(5,0),Transposition(6,2)) );
        BOOST_CHECK( S13.at(2) == s0 );
        BOOST_CHECK( S13.at(3) == s9 );

        T.at(1)=t8;
        T.at(2)=t5;
        T.at(3)=t6;
        BOOST_CHECK(!T.phi(3,1));
    }

    BOOST_FIXTURE_TEST_CASE(test_d_hor_eventuple_parallelcase, TestNorm2Permutations)
    {
        EvenTuple::parallel_case();
        EvenTuple T(6,3);
        EvenTuple boundary1, boundary2, boundary3, boundary4, boundary5;
        T.at(1)=t3;
        T.at(2)=t2;
        T.at(3)=t1;

        boundary1 = T.d_hor_double_complex(1);
        boundary2 = T.d_hor_double_complex(2);
        boundary3 = T.d_hor_double_complex(3);
        boundary4 = T.d_hor_double_complex(4);
        boundary5 = T.d_hor_double_complex(5);
        BOOST_CHECK(boundary1 == EvenTuple());
        BOOST_CHECK(boundary2 != EvenTuple());
        BOOST_CHECK(boundary2.at(1)== Norm2Permutation(Transposition(4,2),Transposition(5,1)));
        BOOST_CHECK(boundary2.at(2)== Norm2Permutation(Transposition(3,2),Transposition(4,1)));
        BOOST_CHECK(boundary2.at(3)== Norm2Permutation(Transposition(2,1),Transposition(4,2)));
        BOOST_CHECK(boundary3 == EvenTuple());
        BOOST_CHECK(boundary4 == EvenTuple());
        BOOST_CHECK(boundary5 != EvenTuple());
        BOOST_CHECK(boundary5.at(1)== Norm2Permutation(Transposition(3,1),Transposition(5,1)));
        BOOST_CHECK(boundary5.at(2)== s2);
        BOOST_CHECK(boundary5.at(3)== s1);

        boundary2 =T.d_hor(2);
        boundary5 =T.d_hor(5);
        BOOST_CHECK(boundary2 != EvenTuple());
        BOOST_CHECK(boundary2.at(1)== Norm2Permutation(Transposition(4,2),Transposition(5,1)));
        BOOST_CHECK(boundary2.at(2)== Norm2Permutation(Transposition(3,2),Transposition(4,1)));
        BOOST_CHECK(boundary2.at(3)== Norm2Permutation(Transposition(2,1),Transposition(4,2)));
        BOOST_CHECK(boundary5 != EvenTuple());
        BOOST_CHECK(boundary5.at(1)== Norm2Permutation(Transposition(3,1),Transposition(5,1)));
        BOOST_CHECK(boundary5.at(2)== s2);
        BOOST_CHECK(boundary5.at(3)== s1 );


        EvenTuple S(6,1), boundary;
        S.at(1) = Norm2Permutation(Transposition(2,1),Transposition(3,2));
        boundary = S.d_hor_double_complex(1);
        BOOST_CHECK(boundary == EvenTuple());
        boundary = S.d_hor(1);
        BOOST_CHECK(boundary == EvenTuple());
    }

    BOOST_FIXTURE_TEST_CASE(test_d_hor_eventuple_radialcase, TestNorm2Permutations)
    {
        EvenTuple::radial_case();
        EvenTuple T(6,3);
        EvenTuple boundary0, boundary1, boundary2, boundary3, boundary4, boundary5, boundary6;
        T.at(1)=t3;
        T.at(2)=t2;
        T.at(3)=t1;

        boundary0 = T.d_hor_double_complex(0);
        boundary1 = T.d_hor_double_complex(1);
        boundary2 = T.d_hor_double_complex(2);
        boundary3 = T.d_hor_double_complex(3);
        boundary4 = T.d_hor_double_complex(4);
        boundary5 = T.d_hor_double_complex(5);
        boundary6 = T.d_hor_double_complex(6);
        BOOST_CHECK(boundary0 != EvenTuple());
        BOOST_CHECK(boundary0.at(1)== Norm2Permutation(Transposition(4,2),Transposition(5,0)));
        BOOST_CHECK(boundary0.at(2)== Norm2Permutation(Transposition(1,0),Transposition(3,2)));
        BOOST_CHECK(boundary0.at(3)== Norm2Permutation(Transposition(1,0),Transposition(2,0)));
        BOOST_CHECK(boundary1 == EvenTuple());
        BOOST_CHECK(boundary2 != EvenTuple());
        BOOST_CHECK(boundary2.at(1)== Norm2Permutation(Transposition(4,2),Transposition(5,1)));
        BOOST_CHECK(boundary2.at(2)== Norm2Permutation(Transposition(3,2),Transposition(4,1)));
        BOOST_CHECK(boundary2.at(3)== Norm2Permutation(Transposition(2,1),Transposition(4,2)));
        BOOST_CHECK(boundary3 == EvenTuple());
        BOOST_CHECK(boundary4 == EvenTuple());
        BOOST_CHECK(boundary5 != EvenTuple());
        BOOST_CHECK(boundary5.at(1)== Norm2Permutation(Transposition(3,1),Transposition(5,1)));
        BOOST_CHECK(boundary5.at(2)== s2);
        BOOST_CHECK(boundary5.at(3)== s1);
        BOOST_CHECK(boundary6 != EvenTuple());
        BOOST_CHECK(boundary6.at(1)== Norm2Permutation(Transposition(1,0),Transposition(5,3)));
        BOOST_CHECK(boundary6.at(2)== s2);
        BOOST_CHECK(boundary6.at(3)== s1 );
    }

    BOOST_FIXTURE_TEST_CASE(test_sigma_h_halves_of_eventuple, TestNorm2Permutations)
    {
        EvenTuple::parallel_case();
        EvenTuple T(6,3);
        T.at(1)=t1;
        T.at(2)=t2;
        T.at(3)=t3;

        Permutation sigma = T.sigma_h_halves();
        BOOST_CHECK(sigma[0]==4 and sigma[4]==3 and sigma[3]==5 and sigma[5]==1
                    and sigma[1]==2 and sigma[2]==6 and sigma[6]==0);

        BOOST_CHECK(T.num_cycles()==1);
        BOOST_CHECK(T.has_correct_num_cycles(0));

        EvenTuple::radial_case();
        BOOST_CHECK(T.has_correct_num_cycles(1));
    }


BOOST_AUTO_TEST_SUITE_END()