// The software kappa is a collection of programs to compute the homology of
// the moduli space of surfaces using the radial model.
// Copyright (C) 2013 - 2018  Felix Boes and Anna Hermann
// 
// This file is part of kappa.
// 
// kappa is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// kappa is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with kappa.  If not, see <http://www.gnu.org/licenses/>.


#ifndef SYM_GRP_TUPLE_HPP
#define SYM_GRP_TUPLE_HPP

#include <cinttypes>
#include <cstdarg>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

#include "permutation.hpp"
#include "permutation_manager.hpp"

/**
 *  A Tuple of Transpositions \f$ ( \tau_{num_entries} \mid \ldots \mid \tau_1 ) \f$.
 *
 *  Note: usually (eg. whenever the SymGrpTuple is a generator in the Ehrenfried complex) we are only interested in the
 *  case num_entries = h, where h is 2*g+m resp. 2*g+m-1 in the parallel resp. radial case. The class can be used for
 *  arbitrary num_entries anyways.
 *
 *  Convention: We want an easy-to-check condition to test if a SymGrpTuple is in a bad state.
 *  We require that bad SymGrpTuples don't contain any transposition.
**/
class SymGrpTuple
{
public:
    friend class HashSymGrpTuple;
    
    /**
     *  Data structure to store the connected components.
     *  The zeroth entry stores the number of connected components.
     *  The components are labled by the natural numbers > 0 and the i-th entry of ConnectedComponents is the component in which i lies.
     */
    typedef std::vector<int32_t> ConnectedComponents;

    /**
     *  Construct a SymGrpTuple with zero entries.
     */
    SymGrpTuple();

    /**
     *  Construct a SymGrpTuple of norm num_entries which has to be filled.
     */
    SymGrpTuple(const size_t num_entries);
    
    /**
     *  Construct a SymGrpTuple of norm num_entries on p symbols which has to be filled.
     */
    SymGrpTuple(const uint32_t p, const size_t num_entries);

    /**
     *  Access the \f$ i \f$-th Transposition of the SymGrpTuple. No bounds checked.
     */
    Transposition& operator[](const size_t n);

    /**
     *  Access the \f$ i \f$-th Transposition of the SymGrpTuple (const version).
     */
    const Transposition& at(const size_t n) const;
    
    /**
     *  Access the j-th Transposition.
     */
    Transposition& at(const size_t j);
    
    /**
     *  @return Returns the number of entries (i.e. the number of transpositions).
     */
    int32_t num_entries() const;
    
    /**
     *  @return Returns true iff both SymGrpTuples are elementwise equal.
     */
    bool operator==(const SymGrpTuple& t) const;
    
    /**
     *  @return Returns false iff both SymGrpTuples are elementwise equal.
     */
    bool operator!=(const SymGrpTuple& t) const;
    
    /**
     *  @return Returns true iff the SymGrpTuple is not marked as degenerate.
     */
    operator bool() const;

    /** 
     *  output stream
     */
    friend std::ostream& operator<< (std::ostream& stream, const SymGrpTuple& tuple);
    
    /**
     *  Tells the number of cycles of the permutation
     *    \f$ \sigma_{num_entries} = \tau_{num_entries} \cdot \ldots \cdot \tau_1 \cdot (1,2,\ldots,p) \f$.
     * \param[in] min_symbol Minimum symbol that may be used in this SymGrpTuple. The default value
     * is 1 since this is the case for parallel cells. Radial cells may also use the
     * symbol 0.
     */
    uint32_t num_cycles() const;

    /*!
     * \brief Determines whether this SymGrpTuple has the correct number of cycles.
     * \param m number of punctures
     * \param radial true iff this SymGrpTuple represents a radial cell
     * \return In the parallel case: true iff the number of cycles equals m + 1;
     *         in the radial case: true iff the number of cycles equals m.
     */
    bool has_correct_num_cycles(const size_t m) const;

    /*!
     * \brief Determines wether this SymGrpTuple is in the product of the cocycle a with something else.
     * \returns In the parallel case: true iff it is the product.
     *          in the radial case: not yet implemented correctly.
     */
    bool is_multiple_of_a() const;

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
    SymGrpTuple Q_term( const std::vector< size_t >& shuffle_slit_conf, const size_t num_shuffle_pos, const size_t fuse_pos ) const;
    
    /**
     *  @returns true iff the SymGrpTuple is fully unstable wrt. mueta
     *  Note: also called 'monotone', as the SymGrpTuple is fully unstable iff the sequence of the heights of the
     *  transpositions is monotone.
     */
    bool fully_unstable() const;
    
    /**
     *  Applies the function \f$ f_i \f$ for \f$ 1 \le i < num_entries() \f$ and returns true iff the norm is preserved thereby.
     */
    bool f(const uint32_t i);
    
    /**
     *  Applies the function \f$ \Phi^q_i \f$ and returns true iff the norm is preserved thereby.
     */
    bool phi( const uint32_t q, const uint32_t i);
    
    /** Applies the i-th horizontal boundary  \f$ \partial_i^{\prime \prime} and the projection on the monotonous cells.
     *  Returns an empty SymGrpTuple if the boundary is degenerate (with respect to the Ehrenfriedcomplex), the boundary
     *  SymGrpTuple otherwise.
     *  \note The parameter i has to fulfill 0 < i < p.
     */
    SymGrpTuple d_hor( const uint8_t i ) const;
    
    /**
     *  Computes the boundary mod multiples of a.
     */
    SymGrpTuple d_hor_reduced( const uint8_t i ) const;

    /** Applies the i-th horizontal boundary  \f$ \partial_i^{\prime \prime} without the projection on the monotonous cells.
     *  Returns an empty SymGrpTuple if the boundary is degenerate (with respect to the double complex), the boundary
     *  SymGrpTuple otherwise.
     *  \note The parameter i has to fulfill 0 < i < p.
     */
    SymGrpTuple d_hor_double_complex( const uint8_t i ) const;
    
    /**
     *  @returns the orientation signs \f$ \varepsilon_0, \ldots, \varepsilon_p \f$.
     */
    std::map< uint8_t, int8_t > orientation_sign() const;
    
    /**
     *  Consider a SymGrpTuple as radial cell.
     */ 
    static void radial_case();
    
    /**
     *  Consider a SymGrpTuple as parallel cell.
     */ 
    static void parallel_case();
    
    /**
     *  @returns true iff we consider a SymGrpTuple as radial cell.
     */ 
    static bool is_radial();
    
    /**
     *  @returns the minimal symbol allowed in a SymGrpTuple, i.e. 1 in the parallel case and 0 in the radial case.
     */ 
    static uint32_t get_min_symbol();
    
    /**
     *  @returns the minimal face index of a possibly non-degenerate cell, i.e. 1 in the parallel case and 0 in the
     *  radial case.
     */ 
    static uint32_t get_min_boundary_offset();
    
    /**
     *  @returns the difference between p and the maximal face index of a possibly non-degenerate cell,
     *  i.e. 1 in the parallel case and 0 in the radial case.
     */ 
    static uint32_t get_max_boundary_offset();
    
    /**
     *  @returns the product of two SymGrpTuples in the sense of the product of two slit domains.
     */
    friend SymGrpTuple operator*( const SymGrpTuple& v_2, const SymGrpTuple& v_1 );
    
    /**
     *  @returns the data representation.
    **/
    const std::vector< Transposition >& get_data_rep() const;
    
    uint32_t p; ///< The maximum of the symbols \f$ min_symbol \le p \f$ to be permuted.
    size_t id;  ///< The index of this SymGrpTuple in the basis of the MonoComplex.

protected:
    /**
     *  @returns the cycle 0 -> 1 -> ... -> p-1 -> p -> 0.
     */
    Permutation long_cycle() const;
    
    /**
     *  @returns the cycle 0 -> p -> p-2 -> ... -> 1 -> 0.
     */
    Permutation long_cycle_inv() const;
    
    /**
     *  @returns the permutation sigma_out = \f$ \sigma_{num_entries} = \tau_{num_entries} \cdots \tau_1 \sigma_0 \f$,
     *  that is the permutation of the outgoing boundary curve of the (radial) slit domain.
     */
    Permutation sigma_out() const;

    /**
     * @returns the inverse permutation of sigma_out.
     *  Note: this is easier to compute than sigma_out.
     */
    Permutation sigma_out_inv() const;
    
    //    DATA MEMBERS
    std::vector< Transposition > rep;   ///< Representation of a SymGrpTuple of transpositions \tau_1, ..., \tau_1.
    static bool radial;                 ///< true iff we consider a SymGrpTuple as radial cell.
    static uint32_t min_symbol;         ///< is set by radial_case and parallel_case accordingly.
    static uint32_t min_boundary_offset;    ///< is set by radial_case and parallel_case accordingly.
    static uint32_t max_boundary_offset;    ///< is set by radial_case and parallel_case accordingly.

    friend class boost::serialization::access;
    template <class Archive> void serialize(Archive &ar, const unsigned int) ///< Implements the serialization of SymGrpTuple.
    {
        ar & p & id & rep;
    }
};

/**
 *  Create a SymGrpTuple consisting of num_entries transpositions \f$\tau_{num_entries} \ldots \tau_1\f$.
 */
SymGrpTuple create_cell( const size_t num_entries, ... );

// output stream
std::ostream& operator<< (std::ostream& stream, const SymGrpTuple& tuple);

/**
 *  @returns the product of two SymGrpTuples in the sense of the product of two slit domains.
 */
SymGrpTuple operator*( const SymGrpTuple& v_2, const SymGrpTuple& v_1 );

/**
 *  In order to save SymGrpTuples in a hash table (e.g. in MonoBasis) we need a function object, that hashes SymGrpTuples.
 */ 
class HashSymGrpTuple
{
public:
    size_t operator()( const SymGrpTuple& ) const;
    
    friend class boost::serialization::access;
    
    template <class Archive>
    void serialize(Archive &, const unsigned int) ///< Implements the serialization.
    {
    }
};

#endif // SYM_GRP_TUPLE_HPP
