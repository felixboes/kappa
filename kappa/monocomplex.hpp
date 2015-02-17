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
#include "misc.hpp"
#include "sessionconfig.hpp"
#include "tuple.hpp"


/**
    The MonoBasis keeps track of the basis elements of a module in a MonoComplex.
**/
struct MonoBasis
{
    MonoBasis();
    
    /// Add a basis element.
    uint32_t add_basis_element (Tuple& t);
    
    /// output stream
    friend std::ostream& operator<< (std::ostream& stream, const MonoBasis& mb);
    
    /// Returns the number of basis elements.
    uint64_t size() const;

    /// Returns the index of the Tuple that is stored in the MonoBasis or -1.
    int64_t id_of( const Tuple& t ) const;
    
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

MonoBasis load_parallel_mono_basis( const uint32_t g, const uint32_t m, const int32_t p );

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
void update_differential(MatrixType &           differential,
                         const size_t           row,
                         const size_t           column,
                         const int32_t          parity,
                         const int8_t           i,
                         const int8_t           or_sign,
                         const SignConvention & sign_conv);

/**
 * This template specialization updates the differential of type MatrixBool according to the contribution
 * of tuple and its boundary.
 */
template<>
void update_differential(MatrixBool &           differential,
                         const size_t           row,
                         const size_t           column,
                         const int32_t          parity,
                         const int8_t           i,
                         const int8_t           or_sign,
                         const SignConvention & sign_conv);

template<>
void update_differential(MatrixBoolCSS &        differential,
                         const size_t           row,
                         const size_t           column,
                         const int32_t          parity,
                         const int8_t           i,
                         const int8_t           or_sign,
                         const SignConvention & sign_conv);
                         
/**
    This MonoComplex represents a chain complex which is generated by the monotone tuples of transpositions. 
    The MonoComplex can either consist of parallel or of radial cells, which is marked by a flag.
    One can generate its bases and differentials.
**/
template< class MatrixComplex >
class MonoComplex
{
public:
    typedef typename MatrixComplex::CoefficientType CoefficientType;
    typedef typename MatrixComplex::MatrixType MatrixType;
    typedef typename MatrixComplex::HomologyType HomologyType;
    typedef typename MatrixComplex::DiagonalizerType DiagonalizerType;
    typedef MonoComplex< MatrixComplex > ThisType;

    MonoComplex( const uint32_t genus, const uint32_t num_punctures, SignConvention sgn, const uint32_t number_working_threads, const uint32_t number_remaining_threads);
    /** Recursive function initializing the basis_complex.
        In the call of gen_bases with the parameters l, p and tuple, we assume that the first l transpositions
        containing symbols start_symbol, ..., p are fixed and append all possible transpositions at position l+1, applying 	the function recursively in an appropriate way.
        If l == h, we don't append another transposition since we have completed a possible basis element. 
        We check whether its number of cycles is appropriate for it to be a basis element, and if this is
        the case, we add it to the basis in degree p.
    **/
    void gen_bases( const uint32_t l, const uint32_t p, const uint32_t start_symbol, Tuple& tuple);
    
    /**
     *  computes the boundary of a given Tuple and saves the result in the differential.
     */ 
    void compute_boundary( Tuple & tuple, const uint32_t p, MatrixType & differential);
    
    /**
     *  Generates the p-th differential.
     */
    void gen_differential( const int32_t p );
    
    /**
     *  Apply base changes.
    **/
    void apply_base_changes();
    
    /**
     *  Diagoanlize current differential.
     */
    HomologyType diagonalize_current_differential( const int32_t p, uint32_t max_rank = 0, const bool print_duration = true );
    
    /**
     *  @returns a reference to the current differential.
     */
    MatrixType&       get_current_differential();
    
    /**
     *  @returns a reference to the current differential.
     */
    const MatrixType& get_current_differential() const;
    
    /**
     *  @return number of rows of the current differential
    **/
    size_t num_rows() const;
    
    /**
     *  @return number of columns of the current differential
    **/
    size_t num_cols() const;
    
    /**
     *  Access the diagonalizer.
    **/
    DiagonalizerType&       get_diagonalizer();
    
    /**
     *  Access the diagonalizer.
    **/
    const DiagonalizerType& get_diagonalizer() const;
    
    /**
     *  erases the current differential.
     */
    void erase_current_differential();
    
    /**
     *  Compute the kernel at the \f$n\f$-th spot and the torsion at the \f$(n-1)\f$-th spot.
    **/
    HomologyType compute_current_kernel_and_torsion( const int32_t n );
    
    /**
     *  print a basis to std::out.
     */
    void show_basis( const int32_t p ) const;
//protected:

    uint32_t g;                 ///< genus
    uint32_t m;                 ///< number of punctures
    uint32_t h;                 ///< h = 2*g+m for the parallel case; and h = 2*g+m-1 for the radial case.
    uint32_t num_threads;       ///< number of threads used to construct the differential

    SignConvention sign_conv;   ///< The sign convention.
    MatrixComplex matrix_complex;                         ///< underlying matrix complex of this MonoComplex
    std::map< int32_t, MonoBasis > bases;        ///< basis_complex[n] is the n-th MonoBasis, i.e. the basis of the n-th module of this MonoComplex. 
};

typedef std::vector<Tuple> MonocomplexWork;

template< class MatrixComplex >
void monocomplex_work(MonoComplex<MatrixComplex> & monocomplex, MonocomplexWork & work, const uint32_t p, typename MatrixComplex::MatrixType & differential);

/**
 * @brief Formula for the sign appearing in d_hor.
 */
int32_t sign(const int32_t          parity,
             const int8_t           i,
             const int8_t           or_sign,
             const SignConvention & sign_conv );

template< class CoefficientT >
class MonoCochainField : public VectorField< CoefficientT >
{
public:    
    typedef CoefficientT CoefficientType;
    typedef MonoCochainField< CoefficientType > ThisType;
    typedef typename VectorField< CoefficientType > :: ThisType VectorType;
    typedef typename VectorField< CoefficientType > :: VectorStorageType VectorStorageType;

    MonoCochainField( uint32_t genus, uint32_t num_punct, uint32_t cohom_deg ) :
        VectorType(), g(genus), m(num_punct), p(cohom_deg)
    {
        basis = load_parallel_mono_basis(g, m, p);
        VectorType::resize( basis.size() );
    }
    
    CoefficientType & operator()( const Tuple& t )
    {
        return VectorType::operator()( basis.id_of(t) );
    }

    const CoefficientType& at( const Tuple& t ) const
    {
        return VectorType::at( basis.id_of(t) );
    }
    
    std::string set_name( const std::string& new_name )
    {
        return name = new_name;
    } 
    
    // grant std::ostream access in order to print matrices to ostreams.
    template< class T >
    friend std::ostream& operator<< ( std::ostream& stream, const MonoCochainField<T> & cochain );
    
protected:
    MonoBasis basis;
    uint32_t g;
    uint32_t m;
    uint32_t p;
    std::string name;
};

template< class CoefficientT >
std::ostream& operator<< ( std::ostream& stream, const MonoCochainField< CoefficientT > & cochain )
{
    return stream
        << "name = " << cochain.name
        << ", g = " << cochain.g
        << ", m = " << cochain.m
        << ", p = " << cochain.p
        << ", representing vector = "
        << static_cast< const typename MonoCochainField< CoefficientT >::VectorType & >(cochain);
}

#include "monocomplex.ipp"

#endif // MONOCOMPLEX_H
