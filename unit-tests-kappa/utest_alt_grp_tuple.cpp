#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "kappa/alt_grp_tuple.hpp"

BOOST_AUTO_TEST_SUITE(test_suit_alt_grp_tuple)

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

    BOOST_FIXTURE_TEST_CASE(test_access_and_check_norm, TestNorm2Permutations)
    {
        AltGrpTuple::parallel_case();
        BOOST_CHECK(!AltGrpTuple::is_radial());

        AltGrpTuple T(6,3);
        BOOST_CHECK(T.num_entries()==3);

        T.at(1)=t1;
        T.at(2)=t2;
        T[3]=t3;
        BOOST_CHECK(T.num_entries()==3);
        BOOST_CHECK(T.has_correct_norm());

        T.at(3)=t4;
        BOOST_CHECK(!T.has_correct_norm());
    }

    BOOST_FIXTURE_TEST_CASE(test_norm2notation, TestNorm2Permutations)
    {
        AltGrpTuple T(6,3);
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

    BOOST_FIXTURE_TEST_CASE(test_fully_unstable_parallelcase, TestNorm2Permutations)
    {
        AltGrpTuple::parallel_case();
        AltGrpTuple S1(6,3), S2(6,3), S3(6,3);

        S1.at(1)=t1;
        S1.at(2)=t5;
        S1.at(3)=t3;
        BOOST_CHECK(S1.fully_unstable());

        S2.at(1)=t1;
        S2.at(2)=t3;
        S2.at(3)=t5;
        BOOST_CHECK(S2.fully_unstable());

        S3.at(1)=t5;
        S3.at(2)=t1;
        S3.at(3)=t3;
        BOOST_CHECK(!S3.fully_unstable());
    }

    BOOST_FIXTURE_TEST_CASE(test_fully_unstable_radialcase, TestNorm2Permutations)
    {
        AltGrpTuple::radial_case();
        AltGrpTuple S1(6,3), S2(6,3), S3(6,3);

        S1.at(1)=t0;
        S1.at(2)=t5;
        S1.at(3)=t6;
        BOOST_CHECK(S1.fully_unstable());

        S2.at(1)=t0;
        S2.at(2)=t6;
        S2.at(3)=t5;
        BOOST_CHECK(S2.fully_unstable());

        S3.at(1)=t5;
        S3.at(2)=t0;
        S3.at(3)=t6;
        BOOST_CHECK(not S3.fully_unstable());
    }

    BOOST_FIXTURE_TEST_CASE(test_f_parallelcase, TestNorm2Permutations)
    {
        AltGrpTuple::parallel_case();
        AltGrpTuple S1(6,3), S2(6,3), S3(8,2), S4(6,2);

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

    BOOST_FIXTURE_TEST_CASE(test_f_radialcase, TestNorm2Permutations)
    {
        AltGrpTuple::radial_case();
        BOOST_CHECK(AltGrpTuple::is_radial());
        AltGrpTuple S1(6,3), S2(6,3), S3(8,2), S4(6,2);
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

    BOOST_FIXTURE_TEST_CASE(test_phi_parallelcase, TestNorm2Permutations)
    {
        AltGrpTuple::parallel_case();
        AltGrpTuple S(6,3), S11(6,3), S21(6,3), S32(6,3), S13(6,3), T(6,3), N(6,3);

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

    BOOST_FIXTURE_TEST_CASE(test_phi_radialcase, TestNorm2Permutations)
    {
        AltGrpTuple::radial_case();
        AltGrpTuple S(6,3), S11(6,3), S21(6,3), S32(6,3), S13(6,3), T(6,3), N(6,3);

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

    BOOST_FIXTURE_TEST_CASE(test_d_hor_parallelcase, TestNorm2Permutations)
    {
        AltGrpTuple::parallel_case();
        AltGrpTuple T(6,3);
        AltGrpTuple boundary1, boundary2, boundary3, boundary4, boundary5;
        T.at(1)=t3;
        T.at(2)=t2;
        T.at(3)=t1;

        boundary1 = T.d_hor_double_complex(1);
        boundary2 = T.d_hor_double_complex(2);
        boundary3 = T.d_hor_double_complex(3);
        boundary4 = T.d_hor_double_complex(4);
        boundary5 = T.d_hor_double_complex(5);
        BOOST_CHECK(boundary1 == AltGrpTuple());
        BOOST_CHECK(boundary2 != AltGrpTuple());
        BOOST_CHECK(boundary2.at(1)== Norm2Permutation(Transposition(4,2),Transposition(5,1)));
        BOOST_CHECK(boundary2.at(2)== Norm2Permutation(Transposition(3,2),Transposition(4,1)));
        BOOST_CHECK(boundary2.at(3)== Norm2Permutation(Transposition(2,1),Transposition(4,2)));
        BOOST_CHECK(boundary3 == AltGrpTuple());
        BOOST_CHECK(boundary4 == AltGrpTuple());
        BOOST_CHECK(boundary5 != AltGrpTuple());
        BOOST_CHECK(boundary5.at(1)== Norm2Permutation(Transposition(3,1),Transposition(5,1)));
        BOOST_CHECK(boundary5.at(2)== s2);
        BOOST_CHECK(boundary5.at(3)== s1);

        boundary2 =T.d_hor(2);
        boundary5 =T.d_hor(5);
        BOOST_CHECK(boundary2 != AltGrpTuple());
        BOOST_CHECK(boundary2.at(1)== Norm2Permutation(Transposition(4,2),Transposition(5,1)));
        BOOST_CHECK(boundary2.at(2)== Norm2Permutation(Transposition(3,2),Transposition(4,1)));
        BOOST_CHECK(boundary2.at(3)== Norm2Permutation(Transposition(2,1),Transposition(4,2)));
        BOOST_CHECK(boundary5 != AltGrpTuple());
        BOOST_CHECK(boundary5.at(1)== Norm2Permutation(Transposition(3,1),Transposition(5,1)));
        BOOST_CHECK(boundary5.at(2)== s2);
        BOOST_CHECK(boundary5.at(3)== s1 );


        AltGrpTuple S(6,1), boundary;
        S.at(1) = Norm2Permutation(Transposition(2,1),Transposition(3,2));
        boundary = S.d_hor_double_complex(1);
        BOOST_CHECK(boundary == AltGrpTuple());
        boundary = S.d_hor(1);
        BOOST_CHECK(boundary == AltGrpTuple());
    }

    BOOST_FIXTURE_TEST_CASE(test_d_hor_radialcase, TestNorm2Permutations)
    {
        AltGrpTuple::radial_case();
        AltGrpTuple T(6,3);
        AltGrpTuple boundary0, boundary1, boundary2, boundary3, boundary4, boundary5, boundary6;
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
        BOOST_CHECK(boundary0 != AltGrpTuple());
        BOOST_CHECK(boundary0.at(1)== Norm2Permutation(Transposition(4,2),Transposition(5,0)));
        BOOST_CHECK(boundary0.at(2)== Norm2Permutation(Transposition(1,0),Transposition(3,2)));
        BOOST_CHECK(boundary0.at(3)== Norm2Permutation(Transposition(1,0),Transposition(2,0)));
        BOOST_CHECK(boundary1 == AltGrpTuple());
        BOOST_CHECK(boundary2 != AltGrpTuple());
        BOOST_CHECK(boundary2.at(1)== Norm2Permutation(Transposition(4,2),Transposition(5,1)));
        BOOST_CHECK(boundary2.at(2)== Norm2Permutation(Transposition(3,2),Transposition(4,1)));
        BOOST_CHECK(boundary2.at(3)== Norm2Permutation(Transposition(2,1),Transposition(4,2)));
        BOOST_CHECK(boundary3 == AltGrpTuple());
        BOOST_CHECK(boundary4 == AltGrpTuple());
        BOOST_CHECK(boundary5 != AltGrpTuple());
        BOOST_CHECK(boundary5.at(1)== Norm2Permutation(Transposition(3,1),Transposition(5,1)));
        BOOST_CHECK(boundary5.at(2)== s2);
        BOOST_CHECK(boundary5.at(3)== s1);
        BOOST_CHECK(boundary6 != AltGrpTuple());
        BOOST_CHECK(boundary6.at(1)== Norm2Permutation(Transposition(1,0),Transposition(5,3)));
        BOOST_CHECK(boundary6.at(2)== s2);
        BOOST_CHECK(boundary6.at(3)== s1 );
    }

    BOOST_FIXTURE_TEST_CASE(test_sigma_out, TestNorm2Permutations)
    {
        AltGrpTuple::parallel_case();
        AltGrpTuple T(6,3);
        T.at(1)=t1;
        T.at(2)=t2;
        T.at(3)=t3;

        Permutation sigma = T.sigma_out();
        BOOST_CHECK(sigma[0]==4 and sigma[4]==3 and sigma[3]==5 and sigma[5]==1
                    and sigma[1]==2 and sigma[2]==6 and sigma[6]==0);

        BOOST_CHECK(T.num_cycles()==1);
        BOOST_CHECK(T.has_correct_num_cycles(0));

        AltGrpTuple::radial_case();
        BOOST_CHECK(T.has_correct_num_cycles(1));
    }

    BOOST_FIXTURE_TEST_CASE(test_create_alt_Grp_tuple, TestNorm2Permutations)
    {
        AltGrpTuple T1 = create_alt_grp_tuple(3, 3,5,6,1, 1,2,3,4, 1,3,2,3);
        AltGrpTuple T2(6,3);
        T2.at(1)=t1;
        T2.at(2)=t2;
        T2.at(3)=t3;
        BOOST_CHECK(T1 == T2);
    }

    BOOST_AUTO_TEST_CASE(test_has_no_fixed_points_except_zero)
    {
        AltGrpTuple T1 = create_alt_grp_tuple(3, 3,5,6,1, 1,2,3,4, 1,3,2,3);
        AltGrpTuple T2 = create_alt_grp_tuple(3, 3,5,6,0, 0,2,3,4, 0,3,2,3);
        BOOST_CHECK(T1.has_no_common_fixed_points_except_zero());
        BOOST_CHECK(not T2.has_no_common_fixed_points_except_zero());
    }


BOOST_AUTO_TEST_SUITE_END()