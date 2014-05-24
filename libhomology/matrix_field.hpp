#ifndef MATRIX_FIELD_HPP
#define MATRIX_FIELD_HPP

// Description:
//
// This header defines a matrix type with coefficients in an (arbitrary) field.
// Such a field must posess operators like 'operator/= ()'.
//
// Moreover this header offers matrices with Q and Zm coefficients, namely MatrixQ and MatrixZm.

#include <boost/dynamic_bitset.hpp>
#include <iostream>
#include <list>
#include <vector>

#include "field_coefficients.hpp"
#include "parallelization.hpp"

/**
 *  This template class defines a matrix type subject to the field coeffiecients 'CoefficientT'.
 *  Our implementation mimes the functionality of ublas::matrix but performs faster.
 *  You can access elements via operator() e.g. MatrixField M(10,200) M(0,0) = CoefficientT(-3);
 */
template < class CoefficientT >
class MatrixField
{
public:
    typedef CoefficientT CoefficientType;
    typedef std::vector<CoefficientT> MatrixStorageType;    ///< This realizes the implementation of the data.
    typedef std::pair< size_t, size_t > MatrixEntryType;
    typedef std::list< MatrixEntryType > DiagonalType;
    
    MatrixField();  ///< Creates a \f$ 0 \times 0\f$ matrix.
    /**
     *  Creates a matrix with num_rows rows and num_cols columns.
     *  The entries are determined by the standard constructor of CoefficientT.
     *  @warning If the standard constructor of CoefficientT does not create a coefficient with value zero (e.g. usind standard int types)
     *  the result differ from your imagination.
     *  In order to get a zero matrix you may use clear();
     */
    MatrixField( size_t number_rows, size_t number_cols );  ///< 
    
    /**
     *  This performs a row operation in the Gauss algorithm.
     *  The entry in (row_1, col) is the given entry that is used to erase the entry in (row_2, col).
     */ 
    void row_operation( size_t row_1, size_t row_2, size_t col );
    
    /**
     *  In order to access elements of the matrix you want to use this function.
     *  @return The function returns a reference to the given entry.
     *  @todo throw an exception if necessary i.e. if (i,j) is not a valid entry.
     */ 
    CoefficientT & operator()( size_t i, size_t j );
    
    /**
     *  In order to keep constness, we go the usual way and implement the function at.
     *  You may want to take a look at the at()-methods of the standard containers like std::vector.
     */
    const CoefficientT& at( size_t i, size_t j ) const;

    /**
     *  As our implementation mimes ublas::matrix we use the same (awkward) method to delete a matrix.
     *  In order to do so call resize(0,0);
     */ 
    void resize (size_t size1, size_t size2, bool);
    
    size_t size1() const;   ///< @returns the number of rows.
    size_t size2() const;   ///< @returns the number of columns.
    
    void clear();   ///< Fills every entry with CoefficientT(0);
    
     /**
     *  grant std::ostream access in order to print coefficients to ostreams like 'std::cout << Zm(44) << std::endl;'
     */
    template< class T >
    friend std::ostream& operator<< ( std::ostream& stream, const MatrixField<T> & matrix );

    DiagonalType diagonal;

private:
    MatrixStorageType data; ///< This realizes the data.
    size_t num_rows;    ///< The number of rows.
    size_t num_cols;    ///< The number of columns.
    
    /**
     *  In order to save Zm coefficients we have to grad boost::serialization::access access.
     */ 
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive &ar, const unsigned int ) ///< Implements the serialization.
    {
        ar & num_rows & num_cols & data;
    }
};

template< class CoefficientT >
std::ostream& operator<< ( std::ostream& stream, const MatrixField<CoefficientT> & matrix );

template class MatrixField<Q>;
template class MatrixField<Zm>;

typedef MatrixField<Q> MatrixQ;     ///< This defines Matrices with \f$\mathbb Q\f$ coefficients.
typedef MatrixField<Zm> MatrixZm;   ///< This defines Matrices with \f$\mathbb Z/ m\mathbb Zf$ coefficients.

////////////////////////////////////////////////////////////////////////////////////////////////////


/**
 *  This template class defines a matrix type subject to the field coeffiecients 'CoefficientT'.
 *  Our implementation mimes the functionality of ublas::matrix but performs faster.
 *  You can access elements via operator() e.g. MatrixField M(10,200) M(0,0) = CoefficientT(-3);
 *  It is used to compute the cluster spectral sequence step by step, therefore it consists of
 *  a main matrix corresponding to the zero-th ifferential in the css and a matrix that will play the
 *  role of the first differential.
 */
template < class CoefficientT >
class MatrixFieldCSS
{
public:
    typedef CoefficientT CoefficientType;
    typedef std::vector<CoefficientT> MatrixStorageType;    ///< This realizes the implementation of the data.
    typedef std::pair< size_t, size_t > MatrixEntryType;
    typedef std::list< MatrixEntryType > DiagonalType;
    
    enum MatrixFieldCSSInitialization
    {
        only_main,
        only_secondary,
        both
    };
    
    MatrixFieldCSS();  ///< Creates a \f$ 0 \times 0\f$ matrix.
    /**
     *  Creates a matrix with num_rows rows and num_cols columns.
     *  The entries are determined by the standard constructor of CoefficientT.
     *  @warning If the standard constructor of CoefficientT does not create a coefficient with value zero (e.g. usind standard int types)
     *  the result differ from your imagination.
     *  In order to get a zero matrix you may use clear();
     */
    MatrixFieldCSS( MatrixFieldCSSInitialization ini, size_t num_rows1, size_t num_cols1, size_t num_rows2 = 0, size_t num_cols2 = 0 );
    
    MatrixFieldCSS( size_t num_rows1, size_t num_cols1 ); ///< Equivalent to MatrixFieldCSS( only_main, num_rows1, num_cols1, 0, 0 );
    
    /**
     *  This performs a row operation in the Gauss algorithm.
     *  The entry in (row_1, col) is the given entry that is used to erase the entry in (row_2, col).
     *  It is applied to both the \f$d^0\f$ and the \f$d^1\f$ part and is so to speak
     *  triggerd by the \f$d^0\f$ matrix.
     */ 
    void row_operation( size_t row_1, size_t row_2, size_t col );
    
    /**
     *  In order to access elements of the matrix you want to use this function.
     *  @return The function returns a reference to the given entry.
     *  @todo throw an exception if necessary i.e. if (i,j) is not a valid entry.
     */ 
    CoefficientT & operator()( size_t i, size_t j );
    
    /**
     *  In order to keep constness, we go the usual way and implement the function at.
     *  You may want to take a look at the at()-methods of the standard containers like std::vector.
     */
    const CoefficientT& at( size_t i, size_t j ) const;
    
    /**
     *  In order to access elements of the secondary matrix you want to use this function.
     *  @return The function returns a reference to the given entry.
     *  @todo throw an exception if necessary i.e. if (i,j) is not a valid entry.
     */ 
    CoefficientT & sec( size_t i, size_t j );
    
    /**
     *  In order to keep constness, we go the usual way and implement the function at.
     *  You may want to take a look at the at()-methods of the standard containers like std::vector.
     */
    const CoefficientT& sec_at( size_t i, size_t j ) const;
    
    /**
     *  As our implementation mimes ublas::matrix we use the same (awkward) method to delete a matrix.
     *  In order to do so call resize(0,0);
     */ 
    void resize (size_t size1, size_t size2, bool);
    
    /**
     *  As our implementation mimes ublas::matrix we use the same (awkward) method to delete a matrix.
     *  In order to do so call sec_resize(0,0);
     */ 
    void sec_resize (size_t size1, size_t size2, bool);
    
    size_t size1() const;   ///< @returns the number of rows.
    size_t size2() const;   ///< @returns the number of columns.
    
    size_t sec_size1() const;   ///< @returns the number of rows of the secondary matrix.
    size_t sec_size2() const;   ///< @returns the number of columns of the secondary matrix.
    
    void clear();   ///< Fills every entry with CoefficientT(0);
    void sec_clear();   ///< Fills every entry of the secondary matrix with CoefficientT(0);
    
     /**
     *  grant std::ostream access in order to print matrix to ostreams.
     */
    template< class T >
    friend std::ostream& operator<< ( std::ostream& stream, const MatrixFieldCSS<T> & matrix );

    DiagonalType diagonal;

private:    
    MatrixStorageType data; ///< This realizes the data.
    MatrixStorageType sec_data; ///< This realizes the data.
    
    size_t num_rows;    ///< The number of rows.
    size_t num_cols;    ///< The number of columns.
    size_t sec_num_rows;    ///< The number of rows.
    size_t sec_num_cols;    ///< The number of columns.
    
    /**
     *  In order to save Zm coefficients we have to grad boost::serialization::access access.
     */ 
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive &ar, const unsigned int ) ///< Implements the serialization.
    {
        ar & num_rows & num_cols & data 
           & sec_num_rows & sec_num_cols & sec_data;
    }
};

template< class CoefficientT >
std::ostream& operator<< ( std::ostream& stream, const MatrixFieldCSS<CoefficientT> & matrix );

template class MatrixFieldCSS<Q>;
template class MatrixFieldCSS<Zm>;

typedef MatrixFieldCSS<Q> MatrixCSSQ;     ///< This defines Matrices with \f$\mathbb Q\f$ coefficients.
typedef MatrixFieldCSS<Zm> MatrixCSSZm;   ///< This defines Matrices with \f$\mathbb Z/ m\mathbb Zf$ coefficients.

////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Class for a matrix with boolean coefficients
 * A MatrixBool has the same structure as the MatrixField,
 * but for a better performance, here each row of the matrix is a dynamic_bitsets,
 * and for thread safety, the total matrix is a vector of its rows.
 */
class MatrixBool
{
public:
    typedef bool CoefficientType;
    typedef std::vector<boost::dynamic_bitset<> > MatrixStorageType;    ///< This realizes the implementation of the data.
    typedef std::pair< size_t, size_t > MatrixEntryType;
    typedef std::list< MatrixEntryType > DiagonalType;

    MatrixBool();  ///< Creates a \f$ 0 \times 0\f$ matrix.
    /**
     *  Creates a matrix with num_rows rows and num_cols columns with entries 0.
     *  In order to get a zero matrix you may use clear();
     */
    MatrixBool( size_t number_rows, size_t number_cols );  ///<

    /**
     *  This performs a row operation in the Gauss algorithm.
     *  The entry in (row_1, col) is the given entry that is used to erase the entry in (row_2, col).
     */
    void row_operation( size_t row_1, size_t row_2, size_t col );

    /**
     *  In order to read elements of the matrix you want to use this function.
     *  @return The function returns a copy of the given entry.
     *  @todo throw an exception if necessary i.e. if (i,j) is not a valid entry.
     */
    bool operator()( size_t i, size_t j );
    /**
     *  In order to keep constness, we go the usual way and implement the function at.
     *  You may want to take a look at the at()-methods of the standard containers like std::vector.
     */
    bool at( size_t i, size_t j ) const;

    /**
     * Adds 1 to the entry (i, j) of the matrix.
     * @note We need this function for writing access to the matrix
     * since the dynamic bitsets does not provide a reasonable reference operator.
     */
    void add_entry( size_t i, size_t j);

    /**
     *  As our implementation mimes ublas::matrix we use the same (awkward) method to delete a matrix.
     *  In order to do so call resize(0,0);
     */
    void resize (size_t size1, size_t size2, bool);

    size_t size1() const;   ///< @returns the number of rows.
    size_t size2() const;   ///< @returns the number of columns.

    void clear();   ///< Fills every entry with 0;

     /**
     *  grant std::ostream access in order to print coefficients to ostreams like 'std::cout << Zm(44) << std::endl;'
     */
    friend std::ostream& operator<< ( std::ostream& stream, const MatrixBool & matrix);

    DiagonalType diagonal;

private:
    MatrixStorageType data; ///< This realizes the data.
    size_t num_rows;    ///< The number of rows.
    size_t num_cols;    ///< The number of columns.

    /**
     *  In order to save coefficients we have to grad boost::serialization::access access.
     */
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive &ar, const unsigned int ) ///< Implements the serialization.
    {
        ar & num_rows & num_cols & data;
    }
};

std::ostream& operator<< ( std::ostream& stream, const MatrixBool & matrix);

#include "matrix_field.ipp"
#endif // MATRIX_FIELD_HPP
