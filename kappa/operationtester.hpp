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

    void load_basis( MonoIndex );
    void load_basis( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p );
    
    void forget_basis( MonoIndex );
    void forget_basis( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p );
    
    void load_differential( MonoIndex );
    void load_differential( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p );
    
    void forget_differential( MonoIndex );
    void forget_differential( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p );
    
    void load_triangular( MonoIndex );
    void load_triangular( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p );
    
    void forget_triangular( MonoIndex );
    void forget_triangular( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p );
    
    void load_diagonal( MonoIndex );
    void load_diagonal( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p );
    
    void forget_diagonal( MonoIndex );
    void forget_diagonal( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p );
    
    /**
     *  @returns true iff the number of entries of v is the number of basis elements.
     */
    bool vector_is_valid( const MonoIndex& monoindex, const VectorType& v ) const;
    
    bool vector_is_cycle( const MonoIndex& monoindex, const VectorType& v ) const;
    
    
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

    void load_basis( CSSIndex );
    void load_basis( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p, int32_t l );
    
    void forget_basis( CSSIndex );
    void forget_basis( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p, int32_t l );
    
    void load_differential( CSSIndex );
    void load_differential( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p, int32_t l_dom, int32_t l_cod );
    
    void forget_differential( CSSIndex );
    void forget_differential( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p, int32_t l_dom, int32_t l_cod );
    
    void load_triangular( CSSIndex );
    void load_triangular( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p, int32_t l_dom, int32_t l_cod );
    
    void forget_triangular( CSSIndex );
    void forget_triangular( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p, int32_t l_dom, int32_t l_cod );
    
    void load_diagonal( CSSIndex );
    void load_diagonal( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p, int32_t l_dom, int32_t l_cod );
    
    void forget_diagonal( CSSIndex );
    void forget_diagonal( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p, int32_t l_dom, int32_t l_cod );
    
    /**
     *  @returns true iff the number of entries of v is the number of basis elements.
     */
    bool vector_is_valid( const CSSIndex& monoindex, const VectorType& v ) const;
    bool vector_is_cycle( const CSSIndex& monoindex, const VectorType& v ) const;
    
protected:
    std::map< CSSIndex, CSSBasis > basis;
    std::map< CSSIndex, MatrixType > diff;
};

#endif // OPERATIONTESTER_HPP
