#ifndef HIGH_CELL_HPP
#define HIGH_CELL_HPP

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

#include "tuple.hpp"

/**
 *  A HighCell is a cell with \f$ q = h, h-1 \f$.
 *
 *  Convention: We want an easy-to-check condition to test if a HighCell is in a bad state.
 *  We require that bad HighCells don't contain any transposition.
**/
class HighCell
{
public:
    friend class HashHighCell;
    
    /**
     *  Data structure to store the connected components.
     *  The zeroth entry stores the number of connected components.
     *  The components are labled by the natural numbers > 0 and the i-th entry of ConnectedComponents is the component in which i lies.
     */
    typedef std::vector<int32_t> ConnectedComponents;

    /**
     *  Construct a HighCell of norm zero.
     */
    HighCell();

    /**
     *  Construct a HighCell of norm h which has to be filled.
     */
    HighCell(const size_t h);
    
    /**
     *  Construct a HighCell of norm h on p symbols which has to be filled.
     */
    HighCell(const uint32_t p, const size_t h);

    /**
     *  Access the \f$ i \f$-th Transposition of the HighCell.
     */
    Transposition& operator[](const size_t n);

    /**
     *  Access the \f$ i \f$-th Transposition of the HighCell (const version).
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
     *  @returns true iff the cell is redundant.
    **/
    bool is_redundant() const;
    
    /**
     *  @returns the redundancy index if any.
    **/
    uint32_t get_redundancy_index() const;
    
    /**
     *  @return Returns true iff both HighCells are elementwise equal.
     */
    bool operator==(const HighCell& t) const;
    
    /**
     *  @return Returns false iff both HighCells are elementwise equal.
     */
    bool operator!=(const HighCell& t) const;
    
    /**
     *  @return Returns true iff the HighCell is not marked as degenerate.
     */
    operator bool() const;

    /** 
     *  output stream
     */
    friend std::ostream& operator<< (std::ostream& stream, const HighCell& highcell);
    
    /**
     *  Tells the number of cycles of the permutation \f$ \tau_h \cdot \ldots \cdot \tau_1 \cdot (1,2,\ldots,p) \f$.
     * \param[in] min_symbol Minimum symbol that may be used in this HighCell. The default value
     * is 1 since this is the case for parallel cells. Radial cells may also use the
     * symbol 0.
     */
    uint32_t num_cycles() const;

    /*!
     * \brief Determines whether this HighCell has the correct number of cycles.
     * \param m number of punctures
     * \param radial true iff this HighCell represents a radial cell
     * \return In the parallel case: true iff the number of cycles equals m + 1;
     *         in the radial case: true iff the number of cycles equals m.
     */
    bool has_correct_num_cycles(const size_t m) const;

    /**
     *  @returns the number connected compontents of the corresponding graph, where \f$ \tau_j \f$ is seen an edge.
     */
    ConnectedComponents connected_components() const;
    
    /**
     *  @returns the number of clusters.
     */
    int32_t num_clusters() const;
        
    /**
     *  @returns true iff the HighCell is monotone.
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
    
    /** Applies the i-th horizontal boundary  \f$ \partial_i^{\prime \prime} without the projection on the monotonous cells.
     *  Returns an empty HighCell if the boundary is degenerate (with respect to the double complex), the boundary HighCell otherwise.
     *  \note The parameter i has to fulfill 0 < i < p.
     */
    HighCell d_hor_double_complex( const uint8_t i ) const;
    
    /**
     *  @returns the ith vertical face. 
     */
    HighCell d_ver( const uint8_t i) const;
    
    /**
     *  @returns the orientation signs \f$ \varepsilon_0, \ldots, \varepsilon_p \f$.
     */
    std::map< uint8_t, int8_t > orientation_sign() const;
    
    /**
     *  @returns the minimal symbol allowed in a HighCell, i.e. 1 in the parallel case and 0 in the radial case.
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
     *  @returns the product of two HighCells in the sense of the product of two slit domains.
     */
    friend HighCell operator*( const HighCell& v_2, const HighCell& v_1 );
    
    uint32_t p; ///< The number of symbols \f$ 1 \le p \f$ to be permuted.
    uint32_t redundancy_index;           ///< This is 0 iff the cell is not redundant. In the other case it is paired with exactly one cell in the top dimension by a destinguished coboundary, whose index is stored.
    size_t id;  ///< The index of this HighCell in the basis of the MonoComplex.
    
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
    std::vector< Transposition > rep;   ///< Representation of a HighCell of transpositions \tau_1, ..., \tau_1.
    static const uint32_t min_symbol;   ///< is set by radial_case and parallel_case accordingly.
    static const uint32_t min_boundary_offset;  ///< is set by radial_case and parallel_case accordingly.
    static const uint32_t max_boundary_offset;  ///< is set by radial_case and parallel_case accordingly.

    friend class boost::serialization::access;
    template <class Archive> void serialize(Archive &ar, const unsigned int) ///< Implements the serialization of HighCell.
    {
        ar & p & id & rep;
    }
};

/**
 *  Create a HighCell consisting of h transpositions \f$\tau_h \ldots \tau_1\f$.
 */
HighCell create_highcell( const size_t h, ... );

// output stream
std::ostream& operator<< (std::ostream& stream, const HighCell& highcell);

/**
 *  @returns the product of two HighCells in the sense of the product of two slit domains.
 */
HighCell operator*( const HighCell& v_2, const HighCell& v_1 );

/**
 *  In order to save HighCells in a hash table (e.g. in MonoBasis) we need a function object, that hashes HighCells.
 */ 
class HashHighCell
{
public:
    size_t operator()( const HighCell& ) const;
    
    friend class boost::serialization::access;
    
    template <class Archive>
    void serialize(Archive &, const unsigned int) ///< Implements the serialization.
    {
    }
};

#endif // HIGH_CELL_HPP
