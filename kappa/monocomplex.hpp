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
#include "monobasis.hpp"
#include "operationtester.hpp"
#include "sessionconfig.hpp"
#include "tuple.hpp"


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
    void compute_boundary_sage( Tuple & tuple, const uint32_t p, SagemathInterface& sage);
    
    /**
     *  Generates the p-th differential.
     */
    void gen_differential( const int32_t p );
    void gen_differential_sage( const int32_t p, SagemathInterface& sage);
    
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
        const auto res = basis.id_of(t);
        if( res == -1 )
        {
            std::cout << "Error: " << t << " is no basis element." << std::endl;
        }
        return VectorType::operator()( res );
    }

    const CoefficientType& at( const Tuple& t ) const
    {
        return VectorType::at( basis.id_of(t) );
    }
    
    std::string set_name( const std::string& new_name )
    {
        return name = new_name;
    } 
    
    void add_kappa_dual( const CoefficientType& c, const Tuple& t )
    {
        VectorType::operator+=( kappa_dual< VectorType >( c, t, basis ) );
    } 
    
    template< class T >
    friend ThisType operator*( const ThisType&, const ThisType& );
    
    uint32_t get_g() const
    {
        return g;
    }
    
    uint32_t get_m() const
    {
        return m;
    }
        
    uint32_t get_p() const
    {
        return p;
    }
    
    std::string get_name() const
    {
        return name;
    }
    
    const MonoBasis& get_basis_reference() const
    {
        return basis;
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
MonoCochainField< CoefficientT > operator*( const MonoCochainField< CoefficientT >&  x, const MonoCochainField< CoefficientT >& y )
{
    typedef MonoCochainField< CoefficientT > CochainType;
    CochainType res( x.get_g() + y.get_g(), x.get_m() + y.get_m(), x.get_p() + y.get_p() );
    
    res.set_name( x.get_name() + y.get_name() );
    
    for( const auto cell_x : x.get_basis_reference().basis )
    {
        const CoefficientT& coeff_x = x.at( cell_x );
        if( coeff_x != CoefficientT(0) )
        {
            for( const auto cell_y : y.get_basis_reference().basis )
            {
                const CoefficientT& coeff_y = y.at( cell_y );
                if( coeff_y != CoefficientT(0) )
                {
                    res( cell_x * cell_y ) = coeff_x * coeff_y;
                }
            }
        }
    }
    
    return res;
}

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
