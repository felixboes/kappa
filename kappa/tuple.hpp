#ifndef TUPLE_HPP
#define TUPLE_HPP

#include <cinttypes>
#include <cstdarg>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

#include <libhomology/homology.hpp>

/**
 *  A Transposition is represented by a pair of uint8_t.
 */
typedef std::pair< uint8_t, uint8_t > Transposition;

/**
 * @brief Class representing a permutation
 *
 * The permutation is stored in a vector of a given size, say s,
 * whereby the permutation acts on a subset of the numbers 0, ..., s-1.
 * For each symbol 0 <= i < s, its entry in the vector tells its image under the permutation.
 * Per default, each value i maps to s, marking that i does not actually belong to the permutation.
 * 
 */
class Permutation
{
public:
    /**
     * Constructs an empty permutation.
     */
    Permutation();

    /**
     * Constructs a permutation of the given size and initializes each entry with the default value size.
     */
    Permutation(const uint8_t size);

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
    
    /**
     * Returns the element the symbol i is mapped to by this Permutation.
     */
    const uint8_t & at(const uint8_t i) const;

    /**
     * Returns the number of elements of this permutation.
     */
    uint8_t size() const;
    
    /**
     *  @brief Determines the decompositions of the permutation into cycles including fix points.
     *  @return A map consisting of all cycles of pi with their smallest element as key
     */
    std::map< uint8_t, Permutation > cycle_decomposition () const;
    
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

protected:
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
public:
    friend class HashTuple;
    
    /**
     *  Data structure to store the connected components.
     *  The zeroth entry stores the number of connected components.
     *  The components are labled by the natural numbers > 0 and the i-th entry of ConnectedComponents is the component in which i lies.
     */
    typedef std::vector<int32_t> ConnectedComponents;

    /**
     *  Construct a Tuple of norm zero.
     */
    Tuple();

    /**
     *  Construct a Tuple of norm h which has to be filled.
     */
    Tuple(const size_t h);
    
    /**
     *  Construct a Tuple of norm h on p symbols which has to be filled.
     */
    Tuple(const uint32_t p, const size_t h);

    /**
     *  Access the \f$ i \f$-th Transposition of the Tuple.
     */
    Transposition& operator[](const size_t n);

    /**
     *  Access the \f$ i \f$-th Transposition of the Tuple (const version).
     */
    const Transposition& at(const size_t n) const;
    
    /**
     *  Access the j-th Transposition.
     */
    Transposition& at(const size_t j);
    
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

    /**
     *  @returns the number connected compontents of the corresponding graph, where \f$ \tau_j \f$ is seen an edge.
     */
    ConnectedComponents connected_components() const;
    
    /**
     *  @returns the number of clusters.
     */
    int32_t num_clusters() const;
    
    /**
     *  @returns a vector encoding the shuffle positions.
     *  These are the slits occuring by traversing from top to bottom on the left hand side of the cell.
     *  The last entry is a zero.
     *  Here is a little example.
     *  The cell \f$ \big(\ (3,1) \mid (2,1) \big) \f$ has four slits enumerated by their occurence in \f$ tau_j \f$ and their height.
     *  More precisely the upper slit corresponding to \f$ \tau_j\f$ is \f$2j\f$ and the lower is \f$ 2j-1 \f$.
     *  Here we tranverse the slits 4, 1, 3, 2.
    **/
    std::vector< size_t > shuffle_positions() const;
    
    /**
     *  @returns an enumeration of the slits.
     *  Here, we count the slits from bottom to top.
     *  The position k means that the slit occures in the k/2-th transposition e.g. the slits of the cell
     *  \f$ \big(\ (3,1) \mid (2,1) \big) /f$ are enumerated as follows.
     *  \f$ 1, 3, 2, 4\f$.
    **/
    std::vector< size_t > slits() const;
    
    /**
     * @brief generates a single term of Q with respect to the given data.
    **/
    Tuple Q_term( const std::vector< size_t >& shuffle_slit_conf, const size_t num_shuffle_pos, const size_t fuse_pos ) const;
    
    /**
     *  @returns true iff the tuple is monotone.
     */
    bool monotone() const;
    
    /**
     *  Applies the function \f$ f_i \f$ for \f$ 1 \le i < h \f$ and returns true iff the norm is preserved thereby.
     */
    bool f(const uint32_t i);
    
    /**
     *  Applies the function \f$ \Phi^q_i \f$ and returns true iff the norm is preserved thereby.
     */
    bool phi( const uint32_t q, const uint32_t i);
    
    /** Applies the i-th horizontal boundary  \f$ \partial_i^{\prime \prime} and the projection on the monotonous cells.
     *  Returns an empty tuple if the boundary is degenerate (with respect to the Ehrenfriedcomplex), the boundary tuple otherwise.
     *  \note The parameter i has to fulfill 0 < i < p.
     */
    Tuple d_hor( const uint8_t i ) const;
    
    /** Applies the i-th horizontal boundary  \f$ \partial_i^{\prime \prime} without the projection on the monotonous cells.
     *  Returns an empty tuple if the boundary is degenerate (with respect to the double complex), the boundary tuple otherwise.
     *  \note The parameter i has to fulfill 0 < i < p.
     */
    Tuple d_hor_double_complex( const uint8_t i ) const;
    
    /**
     *  @returns the orientation signs \f$ \varepsilon_0, \ldots, \varepsilon_p \f$.
     */
    std::map< uint8_t, int8_t > orientation_sign() const;
    
    /**
     *  Consider a Tuple as radial cell.
     */ 
    static void radial_case();
    
    /**
     *  Consider a Tuple as parallel cell.
     */ 
    static void parallel_case();
    
    /**
     *  @returns true iff we consider a Tuple as radial cell.
     */ 
    static bool get_radial();
    
    /**
     *  @returns the minimal symbol allowed in a Tuple, i.e. 1 in the parallel case and 0 in the radial case.
     */ 
    static uint32_t get_min_symbol();
    
    /**
     *  @returns the minimal face index of a possibly non-degenerate cell, i.e. 1 in the parallel case and 0 in the radial case.
     */ 
    static uint32_t get_min_boundary_offset();
    
    /**
     *  @returns the maximal face index of a possibly non-degenerate cell, i.e. p-1 in the parallel case and p in the radial case.
     */ 
    static uint32_t get_max_boundary_offset();
    
    /**
     *  @returns the product of two tuples in the sense of the product of two slit domains.
     */
    friend Tuple operator*( const Tuple& v_2, const Tuple& v_1 );
    
    uint32_t p; ///< The number of symbols \f$ 1 \le p \f$ to be permuted.
    size_t id;  ///< The index of this Tuple in the basis of the MonoComplex.

protected:
    /**
     *  @returns the cycle 1 -> 2 -> ... -> p-1 -> p -> 1.
     */
    Permutation long_cycle() const;
    
    /**
     *  @returns the cycle 1 -> p -> p-2 -> ... -> 2 -> 1.
     */
    Permutation long_cycle_inv() const;
    
    /**
     *  @returns the permutation \f$ \sigma_h = \tau_h \cdots \tau_1 \sigma_0 \f$.
     */
    Permutation sigma_h() const;
    
    //    DATA MEMBERS
    std::vector< Transposition > rep;   ///< Representation of a tuple of transpositions \tau_1, ..., \tau_1.
    static bool radial;                 ///< true iff we consider a Tuple as radial cell.
    static uint32_t min_symbol;         ///< is set by radial_case and parallel_case accordingly.
    static uint32_t min_boundary_offset;    ///< is set by radial_case and parallel_case accordingly.
    static uint32_t max_boundary_offset;    ///< is set by radial_case and parallel_case accordingly.

    friend class boost::serialization::access;
    template <class Archive> void serialize(Archive &ar, const unsigned int) ///< Implements the serialization of Tuple.
    {
        ar & p & id & rep;
    }
};

/**
 *  Create a tuple consisting of h transpositions \f$\tau_h \ldots \tau_1\f$.
 */
Tuple create_cell( const size_t h, ... );

// output stream
std::ostream& operator<< (std::ostream& stream, const Tuple& tuple);

/**
 *  @returns the product of two tuples in the sense of the product of two slit domains.
 */
Tuple operator*( const Tuple& v_2, const Tuple& v_1 );

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
