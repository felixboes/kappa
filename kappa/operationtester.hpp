#ifndef OPERATIONTESTER_HPP
#define OPERATIONTESTER_HPP

#include <string>
#include <tuple>

#include <boost/range/adaptor/reversed.hpp>

#include <libhomology/homology.hpp>

#include "factorial.hpp"
#include "misc.hpp"
#include "monobasis.hpp"
#include "cssbasis.hpp"

template< class MatrixComplex, class VectorT >
class OperationTester
{
public:
    typedef OperationTester< MatrixComplex, VectorT > ThisType;
    typedef typename MatrixComplex::CoefficientType CoefficientType;
    typedef typename MatrixComplex::MatrixType MatrixType;
    typedef typename MatrixType::DiagonalType DiagonalType;
    typedef typename MatrixComplex::HomologyType HomologyType;
    typedef typename MatrixComplex::DiagonalizerType DiagonalizerType;
    typedef VectorT VectorType;
    // radial, g, m, p
    typedef std::tuple< bool, uint32_t, uint32_t, int32_t  > MonoIndex;
    
    OperationTester( std::string coeff_prefix );

    /**
     *  Tries to load a given basis.
     *  @returns true on success.
     */
    bool load_basis( const MonoIndex& idx, bool print_status_messages = true );
    /**
     *  Tries to load a given basis.
     *  @returns true on success.
     */
    bool load_basis( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p, bool print_status_messages = true );
    
    /**
     *  @returns the dimension of a given basis.
     */ 
    size_t dim( const MonoIndex& idx ) const;
    /**
     *  @returns the dimension of a given basis.
     */ 
    size_t dim( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p ) const;
    
    /**
     *  Free given data.
     */
    void forget_basis( const MonoIndex& idx );
    /**
     *  Free given data.
     */
    void forget_basis( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p );
    
    /**
     *  Tries to load a given base change.
     *  @returns true on success.
     */
    bool load_base_changes( const MonoIndex& idx, bool print_status_messages = true );
    /**
     *  Tries to load a given base change.
     *  @returns true on success.
     */
    bool load_base_changes( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p, bool print_status_messages = true );
    
    /**
     *  Free given data.
     */
    void forget_base_changes( const MonoIndex& idx );
    /**
     *  Free given data.
     */
    void forget_base_changes( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p );
    
    /**
     *  Tries to load a given triangular matrix.
     *  @returns true on success.
     */
    bool load_triangular( const MonoIndex& idx, bool print_status_messages = true );
    /**
     *  Tries to load a given triangular matrix.
     *  @returns true on success.
     */
    bool load_triangular( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p, bool print_status_messages = true );
    
    /**
     *  Free given data.
     */
    void forget_triangular( const MonoIndex& idx );
    /**
     *  Free given data.
     */
    void forget_triangular( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p );
    
    /**
     *  Tries to load a given diagonal.
     *  @returns true on success.
     */
    bool load_diagonal( const MonoIndex& idx, bool print_status_messages = true );
    /**
     *  Tries to load a given diagonal.
     *  @returns true on success.
     */
    bool load_diagonal( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p, bool print_status_messages = true );
    
    /**
     *  Free given data.
     */
    void forget_diagonal( const MonoIndex& idx );
    /**
     *  Free given data.
     */
    void forget_diagonal( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p );
    
    /**
     *  Prints status messages corresponding to the cached data.
     */
    void print_cache_status() const;
    
    /**
     *  @returns true iff the number of entries of v is the number of basis elements.
     */
    bool vector_is_valid( const MonoIndex& idx, const VectorType& v ) const;
    
    /**
     *  @returns true iff the chain is a cycle. Tries to load matrix.
     */
    bool vector_is_cycle( const MonoIndex& idx, const VectorType& v );
    
//    /**
//     *  @returns the given homology class of the vector v.
//     *  @warning We assert, that p >= 2 and v is valid cycle.
//     */
//    VectorType vector_homology_class( const MonoIndex& idx, const VectorType& v );
    
    /**
     *  @returns the index of the associated moduli space.
     */
    static MonoIndex  product( const MonoIndex& idx_v, const MonoIndex& idx_w );
    
    /**
     *  @returns the product of two vectors in he associated moduli space.
     */
    VectorType product( const MonoIndex& idx_v, const VectorType& v, const MonoIndex& idx_w, const VectorType& w );
    
//    /**
//     *  @returns Q(v).
//    **/
//    VectorType Q( const MonoIndex& idx_v, const VectorType& v );
    
//    /**
//     *  compute and add Q of a single tuple to a given vector.
//    **/
//    void compute_and_add_Q( const CoefficientType& c, const Tuple& t, const MonoBasis& b, VectorType v );

//protected:
    std::string coefficient_prefix;
    std::map< MonoIndex, MonoBasis > basis;
    std::map< MonoIndex, MatrixType > base_changes;
    std::map< MonoIndex, MatrixType > triangular;
    std::map< MonoIndex, DiagonalType > diagonal;
};

/**
 *  compute and add \f$\kappa^\ast\f$ of a single tuple to a given vector.
**/
template< class VectorT, class CoefficientT >
VectorT kappa_dual( const CoefficientT& c, const Tuple& t, const MonoBasis& b );

template< class VectorT, class CoefficientT >
void compute_and_add_kappa_dual_rec( const CoefficientT& c, const Tuple& t, const MonoBasis& b, VectorT& v, const std::vector<size_t> s, const size_t i );



template< class MatrixComplex, class VectorT >
class OperationTesterCSS
{
public:
    typedef OperationTesterCSS< MatrixComplex, VectorT > ThisType;
    typedef typename MatrixComplex::CoefficientType CoefficientType;
    typedef typename MatrixComplex::MatrixType MatrixType;
    typedef typename MatrixComplex::HomologyType HomologyType;
    typedef typename MatrixComplex::DiagonalizerType DiagonalizerType;
    typedef VectorT VectorType;
    // radial, g, m, p, [(l | l_domain), l_codomain]
    typedef std::tuple< bool, uint32_t, uint32_t, int32_t, int32_t, int32_t  > CSSIndex;
    
    OperationTesterCSS();

    void load_basis( const CSSIndex& idx );
    void load_basis( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p, int32_t l );
    
    size_t dim( const CSSIndex& idx );
    size_t dim( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p );
    
    void forget_basis( const CSSIndex& idx );
    void forget_basis( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p, int32_t l );
    
    bool load_differential( const CSSIndex& idx );
    bool load_differential( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p, int32_t l_dom, int32_t l_cod );
    
    void forget_differential( const CSSIndex& idx );
    void forget_differential( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p, int32_t l_dom, int32_t l_cod );
    
    bool load_triangular( const CSSIndex& idx );
    bool load_triangular( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p, int32_t l_dom, int32_t l_cod );
    
    void forget_triangular( const CSSIndex& idx );
    void forget_triangular( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p, int32_t l_dom, int32_t l_cod );
    
    bool load_diagonal( const CSSIndex& idx );
    bool load_diagonal( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p, int32_t l_dom, int32_t l_cod );
    
    void forget_diagonal( const CSSIndex& idx );
    void forget_diagonal( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p, int32_t l_dom, int32_t l_cod );
    
    /**
     *  @returns true iff the number of entries of v is the number of basis elements.
     */
    bool vector_is_valid( const CSSIndex& monoindex, const VectorType& v ) const;
    bool vector_is_cycle( const CSSIndex& monoindex, const VectorType& v ) const;
    
    /**
     *  @returns the index of the associated moduli space.
     */
    static CSSIndex   product( const CSSIndex& idx_v, const CSSIndex& idx_w );
    
    /**
     *  @returns the vector of the associated moduli space.
     */
    static VectorType product( const CSSIndex& idx_v, const VectorType& v, const CSSIndex& idx_w, const VectorType& w );
    
protected:
    std::map< CSSIndex, CSSBasis > basis;
    std::map< CSSIndex, MatrixType > diff;
};

#endif // OPERATIONTESTER_HPP
