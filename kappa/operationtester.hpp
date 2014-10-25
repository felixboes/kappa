#ifndef OPERATIONTESTER_HPP
#define OPERATIONTESTER_HPP

#include <tuple>

#include "kappa.hpp"

template< class MatrixComplex, class VectorT >
class OperationTester
{
public:
    typedef OperationTester< MatrixComplex > ThisType;
    typedef typename MatrixComplex::CoefficientType CoefficientType;
    typedef typename MatrixComplex::MatrixType MatrixType;
    typedef typename MatrixComplex::HomologyType HomologyType;
    typedef typename MatrixComplex::DiagonalizerType DiagonalizerType;
    typedef VectorT VectorType;
    // radial, g, m, p
    typedef std::tuple< bool, uint32_t, uint32_t, int32_t  > MonoIndex;
    
    OperationTester();

    load_basis( MonoIndex );
    load_basis( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p );
    
    forget_basis( MonoIndex );
    forget_basis( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p );
    
    load_differential( MonoIndex );
    load_differential( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p );
    
    forget_differential( MonoIndex );
    forget_differential( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p );
    
    load_triangular( MonoIndex );
    load_triangular( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p );
    
    forget_triangular( MonoIndex );
    forget_triangular( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p );
    
    load_diagonal( MonoIndex );
    load_diagonal( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p );
    
    forget_diagonal( MonoIndex );
    forget_diagonal( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p );
    
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
    typedef OperationTesterCSS< MatrixComplex > ThisType;
    typedef typename MatrixComplex::CoefficientType CoefficientType;
    typedef typename MatrixComplex::MatrixType MatrixType;
    typedef typename MatrixComplex::HomologyType HomologyType;
    typedef typename MatrixComplex::DiagonalizerType DiagonalizerType;
    typedef VectorT VectorType;
    // radial, g, m, p, [(l | l_domain), l_codomain]
    typedef std::tuple< bool, uint32_t, uint32_t, int32_t, int32_t, int32_t  > CSSIndex;
    
    OperationTesterCSS();

    load_basis( CSSIndex );
    load_basis( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p, int32_t l );
    
    forget_basis( CSSIndex );
    forget_basis( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p, int32_t l );
    
    load_differential( CSSIndex );
    load_differential( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p, int32_t l_dom, int32_t l_cod );
    
    forget_differential( CSSIndex );
    forget_differential( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p, int32_t l_dom, int32_t l_cod );
    
    load_triangular( CSSIndex );
    load_triangular( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p, int32_t l_dom, int32_t l_cod );
    
    forget_triangular( CSSIndex );
    forget_triangular( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p, int32_t l_dom, int32_t l_cod );
    
    load_diagonal( CSSIndex );
    load_diagonal( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p, int32_t l_dom, int32_t l_cod );
    
    forget_diagonal( CSSIndex );
    forget_diagonal( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p, int32_t l_dom, int32_t l_cod );
    
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
