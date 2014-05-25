#ifndef MONOCOMPLEX_H
#define MONOCOMPLEX_H

#include <chrono>
#include <stdint.h>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <unordered_set>

#include <libhomology/homology.hpp>

#include "factorial.hpp"
#include "serialization.hpp"
#include "sessionconfig.hpp"
#include "tuple.hpp"


/**
    The MonoBasis keeps track of the basis elements of a module in a MonoComplex.
**/
struct MonoBasis
{
    /// Add a basis element.
    uint32_t add_basis_element (Tuple& t);
    
    /// output stream
    friend std::ostream& operator<< (std::ostream& stream, const MonoBasis& mb);
    
    /// Returns the number of basis elements.
    uint64_t size() const;

    /// Returns the index of the Tuple that is stored in the MonoBasis or -1.
    int64_t id_of( Tuple& t ) const;
    
    /// Stores the orderd basis.
    std::unordered_set< Tuple, HashTuple > basis;
    
    friend class boost::serialization::access;
    
    /// @warning The serialization library from boost does not yet support unorderd maps (we use boost in the version 1.49). Therefore we must provide a workaround.
    /// boost::serialization method that we use to save a MonoBasis to file.
    template<class Archive>
    void save(Archive & ar, const unsigned int) const
    {
        // In order to load an unorderd_set we need to know the exact number of elemets that are stored.
        size_t size = basis.size();
        ar & size;
        for( const auto& it : basis )
        {
            ar & it;
        }
    }
    
    /// boost::serialization method that we use to load a MonoBasis from file.
    template<class Archive>
    void load(Archive & ar, const unsigned int)
    {
        size_t size;
        Tuple t;
        
        ar & size;
        for( size_t i = 0; i < size; ++i )
        {
            ar & t;
            basis.insert(t);
        }
    }
    
    // This is required as saving and loading are different methods.
    BOOST_SERIALIZATION_SPLIT_MEMBER()
};

std::ostream& operator<< (std::ostream& stream, const MonoBasis& basis);

/**
 * This function is supposed to update the differential according to the contribution
 * of tuple and its boundary to improve readability of the method gen_differential.
 *
 * @warning We don't give a general implementation for this template function since the
 * two specializations we do implement below use specific functions not both of
 * them support.
 *
 * @param parity Parity due to kappa.
 * @param i Number of the differential.
 * @param or_sign Orientation sign.
 */
template <class MatrixType>
void update_differential(MatrixType & differential,
                         Tuple &           tuple,
                         Tuple &           boundary,
                         int32_t           parity,
                         int8_t            i,
                         int8_t            or_sign,
                         SignConvention &  sign_conv);

/**
 * This template specialization updates the differential of type MatrixBool according to the contribution
 * of tuple and its boundary.
 */
template<>
void update_differential(MatrixBool &     differential,
                         Tuple &          tuple,
                         Tuple &          boundary,
                         int32_t          parity,
                         int8_t           i,
                         int8_t           or_sign,
                         SignConvention & sign_conv);

/**
    This MonoComplex represents a chain complex which is generated by the monotone tuples of transpositions. 
**/
template< class MatrixComplex >
class MonoComplex
{
public:
    typedef typename MatrixComplex::CoefficientType CoefficientType;
    typedef typename MatrixComplex::MatrixType MatrixType;
    typedef typename MatrixComplex::HomologyType HomologyType;

    MonoComplex(uint32_t genus, uint32_t num_punctures, SignConvention sgn, uint32_t number_threads, bool radial = false);
    /** Recursive function initializing the basis_complex.
        In the call of gen_bases with the parameters l, p and tuple, we assume that the first l transpositions
        containing symbols min_symbol, ..., p are fixed and append all possible transpositions at position l+1, applying 	the function recursively in an appropriate way.
        If l == h, we don't append another transposition since we have completed a possible basis element. 
        We check whether its number of cycles is appropriate for it to be a basis element, and if this is
        the case, we add it to the basis in degree p.
    **/
    void gen_bases(uint32_t l, uint32_t p, uint32_t min_symbol, Tuple& tuple);
    void compute_boundary(Tuple & tuple, uint32_t p, MatrixType & differential);
    /** Generates the p-th differential.
     * @warning We assume the p-th differential to exist and to be filled with zeroes before the call.
     */
    void gen_differential( int32_t p );
    void gen_differential_naive(int32_t p);                 ///< generate the p-th (the naive way)
    void erase_differential();                              ///< erases the current differential.
    
    void pi_del_phi_naive(const Tuple& it, std::vector<int32_t> & s);
    
    void show_basis( int32_t p ) const;     ///< print a basis to std::out
    void show_differential( int32_t p ) const;  ///< print a differential to std::out
    void show_differential_naive( int32_t p ) const;    ///< print a differential to std::out (if the differential was generated the naive way)
    //std::string show_bases() const;
//private:

    uint32_t g;     ///< genus
    uint32_t m;     ///< number of punctures
    uint32_t h;     ///< h = 2*g+m
    uint32_t num_threads; ///< number of threads used to construct the differential
    bool     radial; ///< true iff this tuple stores a radial cell

    SignConvention sign_conv;  ///< The sign convention.
    MatrixComplex matrix_complex;                         ///< underlying matrix complex of this MonoComplex
    MatrixComplex matrix_complex_naive;                   ///< underlying matrix complex of this MonoComplex (genereted the naive way)
    std::map< int32_t, MonoBasis > basis_complex;        ///< basis_complex[n] is the n-th MonoBasis, i.e. the basis of the n-th module of this MonoComplex. 
};

/**
 * @brief Formula for the sign appearing in d_hor.
 */
int8_t sign(int32_t          parity,
            int8_t           i,
            int8_t           or_sign,
            SignConvention & sign_conv );

#include "monocomplex.ipp"

#endif // MONOCOMPLEX_H
