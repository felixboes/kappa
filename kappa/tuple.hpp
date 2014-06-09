#ifndef TUPLE_HPP
#define TUPLE_HPP

#include <cinttypes>
#include <functional>
#include <iostream>
#include <map>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

#include "serialization.hpp"

/**
 *  A Transposition is represented by a pair of uint8_t.
 */
typedef std::pair< uint8_t, uint8_t > Transposition;

/**
 * @brief Class representing a permutation
 *
 * We assume that the permutation acts on the symbols 1, ..., p.
 * The permutation is stored in the vector data, where for each symbol i in 0, ... p,
 * 0 maps to the symbol data[i]. Thereby, 0 is NOT really part of the permutation!
 * Thus the 0th index of the vector is not used, and we can use data[i] = 0 to
 * represent that the symbol i does not belong to the permutation.
 */
class Permutation
{
public:
    /**
     * Constructs an empty permutation.
     */
    Permutation();

    /**
     * Constructs a permutation of the given size with the given
     * initialization value in each entry.
     */
    Permutation(const uint8_t size, const uint8_t init);

    /**
     * Copy constructor copying the data vector from other.
     */
    Permutation(const Permutation & other);

    /**
     * Returns the data vector of this Permutation.
     */
    std::vector<uint8_t> operator()() const;

    /**
     * Returns the element the symbol i is mapped to by this Permutation.
     */
    uint8_t & operator[](const uint8_t i);
    uint8_t const & at(const uint8_t i) const;

    /**
     * Returns the number of elements of this permutation.
     */
    uint8_t size() const;

    /**
     * Returns true iff the symbol i is contained in this Permutation.
     */
    bool is_contained(const uint8_t i) const;

    /**
     * Returns true iff the symbol i is a fix point of this Permutation.
     */
    bool is_fix_point(const uint8_t i) const;

    /** 
     *  output stream
     */
    friend std::ostream& operator<< (std::ostream& stream, const Permutation& permutation);

private:
    std::vector<uint8_t> data; ///< stores the Permutation
    operator size_t() = delete;
};

/**
 *  A Tuple of Transpositions \f$ ( \tau_h \mid \ldots \mid \tau_1 ) \f$.
 *
 *  Convention: We want an easy-to-check condition to test if a Tuple is in a bad state.
 *  We require that bad Tuples don't contain any transposition.
**/
class Tuple
{
    friend class HashTuple;
    
public:

    /**
     *  Data structure to store the connected components.
     *  The zeroth entry stores the number of connected components.
     *  The components are labled by the natural numbers > 0 and the i-th entry of ConnectedComponents is the component in which i lies.
     */
    
    typedef std::vector<int32_t> ConnectedComponents;

    Tuple();

    /**
     *  Construct a Tuple of norm h which has to be filled.
     */
    Tuple(const size_t h);
    
    Tuple(const uint32_t symbols, const size_t h);

    /**
     *  Access the \f$ i \f$-th Transposition of the Tuple.
     */
    Transposition& operator[](const size_t n);

    /**
     *  @return Returns the norm (i.e. the number of transpositions).
     */
    int32_t norm() const;
    
    /**
     *  @return Returns true iff both Tuples are elementwise equal.
     */
    bool operator==(const Tuple& t) const;
    
    /**
     *  @return Returns false iff both Tuples are elementwise equal.
     */
    bool operator!=(const Tuple& t) const;
    
    /**
     *  @return Returns true iff the tuple is not marked as degenerate.
     */
    operator bool() const;

    /** 
     *  output stream
     */
    friend std::ostream& operator<< (std::ostream& stream, const Tuple& tuple);
    
    /**
     *  Tells the number of cycles of the permutation \f$ \tau_h \cdot \ldots \cdot \tau_1 \cdot (1,2,\ldots,p) \f$.
     * \param[in] min_symbol Minimum symbol that may be used in this tuple. The default value
     * is 1 since this is the case for parallel cells. Radial cells may also use the
     * symbol 0.
     */
    uint32_t num_cycles() const;

    /*!
     * \brief Determines whether this Tuple has the correct number of cycles.
     * \param m number of punctures
     * \param radial true iff this Tuple represents a radial cell
     * \return In the parallel case: true iff the number of cycles equals m + 1;
     *         in the radial case: true iff the number of cycles equals m.
     */
    bool has_correct_num_cycles(size_t m) const;

    ConnectedComponents connected_components() const; ///< @returns the number connected compontents of the corresponding graph, where \f$ \tau_j \f$ is seen an edge.
    int32_t num_clusters() const;
    
    bool monotone();                        ///< Returns true iff the tuple is monotone.
    bool f(const uint32_t i);                     ///< Applies the function \f$ f_i \f$ fuer \f$ 1 \le i < h \f$ and returns true iff the norm is preserved thereby.
    bool phi( const uint32_t q, const uint32_t i);      ///< Applies the function \f$ \Phi^q_i \f$ and returns true iff the norm is preserved thereby..
    /** Applies the i-th horizontal boundary  \f$ \partial_i^{\prime \prime} and the projection on the monotonous cells.
     *  Returns an empty tuple if the boundary is degenerate, the boundary tuple otherwise.
     *  \note The parameter i has to fulfill 0 < i < p.
     */
    Tuple d_hor( const uint8_t i ) const;
    std::map< uint8_t, int8_t > orientation_sign() const;

    uint32_t p;  ///< The number of symbols \f$ 1 \le p \f$ to be permuted.
    size_t id; ///< The index of this Tuple in the basis of the MonoComplex.
    
    static void radial_case();
    static void parallel_case();
    static bool get_radial();
    static uint32_t get_min_symbol();
    static uint32_t get_min_boundary_offset();
    static uint32_t get_max_boundary_offset();
private:
    Transposition& at(const size_t q);                     ///< Access the q-th Transposition.
    Transposition const & at(const size_t q) const;
    Permutation long_cycle() const;                  ///< Returns the cycle 1 -> 2 -> ... -> p-1 -> p -> 1.
    Permutation long_cycle_inv() const;              ///< Returns the cycle 1 -> p -> p-2 -> ... -> 2 -> 1.
    Permutation sigma_h() const;
    std::map< uint8_t, Permutation > cycle_decomposition ( const Permutation & sigma ) const;
    
    //    DATA MEMBERS
    std::vector< Transposition > rep;                ///< Representation of a tuple of transpositions \tau_1, ..., \tau_1.
    static bool radial;
    static uint32_t min_symbol;
    static uint32_t min_boundary_offset;
    static uint32_t max_boundary_offset;

    friend class boost::serialization::access;
    template <class Archive> void serialize(Archive &ar, const unsigned int) ///< Implements the serialization of Tuple.
    {
        ar & p & id & rep;
    }
};

/// output stream
std::ostream& operator<< (std::ostream& stream, const Tuple& tuple);

/**
 *  In order to save Tuples in a hash table (e.g. in MonoBasis) we need a function object, that hashes Tuples.
 */ 
class HashTuple
{
public:
    size_t operator()( const Tuple& ) const;
    
    friend class boost::serialization::access;
    
    template <class Archive>
    void serialize(Archive &, const unsigned int) ///< Implements the serialization.
    {
    }
};

#endif // TUPLE_HPP
