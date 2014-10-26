#ifndef OPERATIONTESTER_HPP
#define OPERATIONTESTER_HPP

#include <tuple>

#include "kappa.hpp"

template< class MatrixComplex, class VectorT >
class OperationTester
{
public:
    typedef OperationTester< MatrixComplex, VectorT > ThisType;
    typedef typename MatrixComplex::CoefficientType CoefficientType;
    typedef typename MatrixComplex::MatrixType MatrixType;
    typedef typename MatrixComplex::HomologyType HomologyType;
    typedef typename MatrixComplex::DiagonalizerType DiagonalizerType;
    typedef VectorT VectorType;
    // radial, g, m, p
    typedef std::tuple< bool, uint32_t, uint32_t, int32_t  > MonoIndex;
    
    OperationTester();

    bool load_basis( const MonoIndex& idx, bool print_status_messages = true );
    bool load_basis( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p, bool print_status_messages = true );
    
    size_t dim( const MonoIndex& idx ) const;
    size_t dim( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p ) const;
    
    void forget_basis( const MonoIndex& idx );
    void forget_basis( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p );
    
    void load_differential( const MonoIndex& idx );
    void load_differential( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p );
    
    void forget_differential( const MonoIndex& idx );
    void forget_differential( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p );
    
    void load_triangular( const MonoIndex& idx );
    void load_triangular( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p );
    
    void forget_triangular( const MonoIndex& idx );
    void forget_triangular( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p );
    
    void load_diagonal( const MonoIndex& idx );
    void load_diagonal( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p );
    
    void forget_diagonal( const MonoIndex& idx );
    void forget_diagonal( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p );
    
    /**
     *  @returns true iff the number of entries of v is the number of basis elements.
     */
    bool vector_is_valid( const MonoIndex& monoindex, const VectorType& v ) const;
    bool vector_is_cycle( const MonoIndex& monoindex, const VectorType& v ) const;
    
    /**
     *  @returns the index of the associated moduli space.
     */
    static MonoIndex  product( const MonoIndex& idx_v, const MonoIndex& idx_w );
    
    /**
     *  @returns the vector of the associated moduli space.
     */
    static VectorType product( const MonoIndex& idx_v, const VectorType& v, const MonoIndex& idx_w, const VectorType& w );
    
    
protected:
    std::map< MonoIndex, MonoBasis > basis;
    std::map< MonoIndex, MatrixType > diff;
};

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
    
    void load_differential( const CSSIndex& idx );
    void load_differential( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p, int32_t l_dom, int32_t l_cod );
    
    void forget_differential( const CSSIndex& idx );
    void forget_differential( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p, int32_t l_dom, int32_t l_cod );
    
    void load_triangular( const CSSIndex& idx );
    void load_triangular( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p, int32_t l_dom, int32_t l_cod );
    
    void forget_triangular( const CSSIndex& idx );
    void forget_triangular( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p, int32_t l_dom, int32_t l_cod );
    
    void load_diagonal( const CSSIndex& idx );
    void load_diagonal( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p, int32_t l_dom, int32_t l_cod );
    
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

typedef OperationTester< ChainComplexQ,  VectorQ  > OperationTesterQ;
typedef OperationTester< ChainComplexZm, VectorZm > OperationTesterZm;
typedef OperationTesterCSS< ChainComplexQCSS,  VectorQ  > OperationTesterQCSS;
typedef OperationTesterCSS< ChainComplexZmCSS, VectorZm > OperationTesterZmCSS;

#include "operationtester.ipp"

#endif // OPERATIONTESTER_HPP
