#ifndef CSS_H
#define CSS_H

#include <algorithm>
#include <chrono>
#include <stdint.h>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdlib.h>
#include <unordered_set>
#include <libhomology/homology.hpp>

#include "factorial.hpp"
#include "monocomplex.hpp"
#include "serialization.hpp"
#include "sessionconfig.hpp"
#include "tuple.hpp"


/**
    The CSSBasis keeps track of the basis elements of a module in the cluster spectral sequence.
**/
struct CSSBasis
{
    /// Stores the orderd basis where every cell has the same number l of cluster.
    typedef std::unordered_set< Tuple, HashTuple > LBasisType;
    
    /// Stores the orderd basis, sorted by cluster sizes.
    typedef std::map< int32_t , LBasisType > BasisType;
    
    /// Add a basis element.
    int32_t add_basis_element ( Tuple& t );
    
    /// output stream
    friend std::ostream& operator<< (std::ostream& stream, const CSSBasis& cssb);
    
    /// Returns the number of basis elements that have cluster size l.
    int32_t size( const int32_t l ) const;
    
    /// Returns the number of basis elements.
    int32_t total_size() const;

    /// Returns the relative (i.e. the cluster-) index of the Tuple that is stored in the CSSBasis or -1.
    int32_t id_of( Tuple& t ) const;
    
    /// Returns the total index of the Tuple that is stored in the CSSBasis or -1.
    int32_t total_id_of( Tuple& t ) const;
    
    BasisType basis;
    
    friend class boost::serialization::access;
    
    /// @warning The serialization library from boost does not yet support unorderd maps (we use boost in the version 1.49). Therefore we must provide a workaround.
    /// boost::serialization method that we use to save a CSSBasis to file.
    template<class Archive>
    void save(Archive & ar, const unsigned int) const
    {
        // In order to load an unorderd_set we need to know the exact number of elemets that are stored.
        size_t num_cluster_sizes = basis.size(); // this is the number of cluster sizes that occur.
        ar & num_cluster_sizes;
        
        for( const auto& it : basis )
        {
            ar & it.first; // This is a specific l.
            for( const auto& b : it.second )
            {
                ar & b; // This is a basis element with exactly l clusters.
            }
        }
    }
    
    /// boost::serialization method that we use to load a CSSBasis from file.
    template<class Archive>
    void load(Archive & ar, const unsigned int)
    {
        size_t num_cluster_sizes;
        
        ar & num_cluster_sizes;
        for( size_t i = 0; i < num_cluster_sizes; ++i )
        {
            size_t l;
            ar & l;
            Tuple t;
            for( size_t j = 0; j < l; ++l )
            {
                ar & t;
                basis[l].insert(t);
            }
        }
    }
    
    // This is required as saving and loading are different methods.
    BOOST_SERIALIZATION_SPLIT_MEMBER()
};

std::ostream& operator<< (std::ostream& stream, const CSSBasis& basis);

/**
This ClusterSpectralSequence represents E_1 term of the cluster spectral sequence associated with the Ehrenfried complex 
which is generated by the monotone tuples of transpositions. 
**/

template< class MatrixComplex >
class ClusterSpectralSequence
{
public:
    typedef typename MatrixComplex::CoefficientType CoefficientType;
    typedef typename MatrixComplex::MatrixType MatrixType;
    typedef typename MatrixComplex::DiagonalizerType DiagonalizerType;
    typedef typename MatrixComplex::HomologyType HomologyType;
    typedef std::map< int32_t, HomologyType > CSSHomologyType;
    
    ClusterSpectralSequence(
            const uint32_t          genus,
            const uint32_t          num_punctures,
            const SignConvention    sgn,
            const uint32_t          number_working_threads,
            const uint32_t          number_remaining_threads );
    /** Recursive function initializing the basis_complex.
        In the call of gen_bases with the parameters s, p and tuple, we assume that the first s transpositions
        containing symbols 1, ..., p are fixed and append all possible transpositions at position s+1, applying
        the function recursively in an appropriate way.
        If s == h, we don't append another transposition since we have completed a possible basis element. 
        We check whether its number of cycles is appropriate for it to be a basis element, and if this is
        the case, we add it to the basis in degree p. Thereby, basis elements are sorted according to the number
        of clusters.
    **/
    void gen_bases( const uint32_t l, const uint32_t p, const uint32_t start_symbol, Tuple& tuple);
    void gen_d0( int32_t p, int32_t l );
    void gen_d0_boundary(const Tuple & tuple,
                         const int32_t p,
                         const int32_t l,
                         typename MatrixComplex::MatrixType & differential);
    
    void gen_d1_stage_1( const int32_t p, const int32_t l );
    MatrixType gen_d1_row( const int32_t, const int32_t l, const Tuple& basis_element );
    void gen_d1_apply_operations( MatrixType& row );
    void prepare_d1_diag();
    void erase_d0();
    void erase_d1();
    
    void show_basis( const int32_t p ) const;         ///< print a basis to std::out
//private:

    uint32_t g;     ///< genus
    uint32_t m;     ///< number of punctures
    uint32_t h;     ///< h = 2*g+m
    uint32_t num_threads;  ///< number of threads used to construct the differential
    
    SignConvention sign_conv;                    ///< The sign convention.
    MatrixComplex diff_complex;                  ///< Due to RAM limitations, we are working with at most two matrices at a time. Therefore we do not model the whole spectral sequence.
    std::map< int32_t, CSSBasis > basis_complex; ///< basis_complex[n] is the n-th CSSBasis.
    
};

typedef std::vector<Tuple> CSSWork;

template< class MatrixComplex >
void css_work_0(ClusterSpectralSequence<MatrixComplex> & css,
              CSSWork & work,
              const int32_t p,
              const int32_t l,
              typename MatrixComplex::MatrixType & differential
              );

template< class MatrixComplex >
void css_work_1(ClusterSpectralSequence<MatrixComplex> & css,
              CSSWork & work,
              const int32_t p,
              const int32_t l,
              typename MatrixComplex::MatrixType & differential,
              const size_t num_cols,
              const std::vector< size_t >& offset
              );

// Template specializations for bool-matrices.
template<>
void ClusterSpectralSequence< ChainComplexBoolCSS >::gen_d1_apply_operations( MatrixType& row );

template<>
void css_work_1(ClusterSpectralSequence< ChainComplexBoolCSS >& css,
              CSSWork & work,
              const int32_t p,
              const int32_t l,
              ChainComplexBoolCSS::MatrixType& differential,
              const size_t num_cols,
              const std::vector< size_t >& offset
              );

typedef ClusterSpectralSequence<ChainComplexQCSS> ClusterSpectralSequenceQ;
typedef ClusterSpectralSequence<ChainComplexZmCSS> ClusterSpectralSequenceZm;
typedef ClusterSpectralSequence<ChainComplexBoolCSS> ClusterSpectralSequenceBool;
typedef ClusterSpectralSequence<ChainComplexZStorageOnly> ClusterSpectralSequenceZStorageOnly;

#include "css.ipp"

#endif // CSS_HP
