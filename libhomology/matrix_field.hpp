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
    MatrixField( const size_t number_rows, const size_t number_cols );  ///< 
    
    /**
     *  This performs a row operation in the Gauss algorithm.
     *  The entry in (row_1, col) is the given entry that is used to erase the entry in (row_2, col).
     */ 
    void row_operation( const size_t row_1, const size_t row_2, const size_t col );
    
    /**
     *  In order to access elements of the matrix you want to use this function.
     *  @return The function returns a reference to the given entry.
     *  @todo throw an exception if necessary i.e. if (i,j) is not a valid entry.
     */ 
    CoefficientT & operator()( const size_t i, const size_t j );
    
    /**
     *  In order to keep constness, we go the usual way and implement the function at.
     *  You may want to take a look at the at()-methods of the standard containers like std::vector.
     */
    const CoefficientT& at( const size_t i, const size_t j ) const;

    /**
     *  As our implementation mimes ublas::matrix we use the same (awkward) method to delete a matrix.
     *  In order to do so call resize(0,0);
     */ 
    void resize ( const size_t size1, const size_t size2, const bool );
    
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
    typedef std::vector< CoefficientT > MatrixStorageType;    ///< This realizes the implementation of the data.
    typedef std::pair< size_t, size_t > MatrixEntryType;
    typedef std::list< MatrixEntryType > DiagonalType;
    typedef MatrixFieldCSS< CoefficientType > ThisType;
    
    enum MatrixFieldCSSInitialization
    {
        only_main,
        only_secondary,
        both
    };
    
    enum OperationType {
        main_and_secondary,
        secondary
    };
    
    MatrixFieldCSS();  ///< Creates a \f$ 0 \times 0\f$ matrix.
    /**
     *  Creates a matrix with num_rows rows and num_cols columns.
     *  The entries are determined by the standard constructor of CoefficientT.
     *  @warning If the standard constructor of CoefficientT does not create a coefficient with value zero (e.g. usind standard int types)
     *  the result differ from your imagination.
     *  In order to get a zero matrix you may use clear();
     */
    MatrixFieldCSS( const MatrixFieldCSSInitialization ini, const size_t num_rows1, const size_t num_cols1, const size_t num_rows2 = 0, const size_t num_cols2 = 0 );
    
    MatrixFieldCSS( const size_t num_rows1, const size_t num_cols1 ); ///< Equivalent to MatrixFieldCSS( only_main, num_rows1, num_cols1, 0, 0 );
    
    /**
     *  This performs a row operation in the Gauss algorithm.
     *  The entry in (row_1, col) is the given entry that is used to erase the entry in (row_2, col).
     *  It is applied to both the \f$d^0\f$ and the \f$d^1\f$ part and is so to speak
     *  triggerd by the \f$d^0\f$ matrix.
     */ 
    void row_operation( const size_t row_1, const size_t row_2, const size_t col );
    void row_operation_main_and_secondary( const size_t row_1, const size_t row_2, const size_t col );
    void row_operation_secondary( const size_t row_1, const size_t row_2, const size_t col );
    void define_operations( const OperationType );
    
    /**
     *  In order to access elements of the matrix you want to use this function.
     *  @return The function returns a reference to the given entry.
     *  @todo throw an exception if necessary i.e. if (i,j) is not a valid entry.
     */ 
    CoefficientType & operator()( const size_t i, const size_t j ); 
    CoefficientType & main_op( const size_t i, const size_t j );
    CoefficientType & sec_op( const size_t i, const size_t j );
    
    /**
     *  In order to keep constness, we go the usual way and implement the function at.
     *  You may want to take a look at the at()-methods of the standard containers like std::vector.
     */
    const CoefficientType& at( const size_t i, const size_t j ) const;
    const CoefficientType& main_at ( const size_t i, const size_t j) const;
    const CoefficientType& sec_at( const size_t i, const size_t j ) const;
    
    /**
     *  As our implementation mimes ublas::matrix we use the same (awkward) method to delete a matrix.
     *  In order to do so call resize(0,0);
     */ 
    void resize ( const size_t size1, const size_t size2, const bool );
    
    /**
     *  As our implementation mimes ublas::matrix we use the same (awkward) method to delete a matrix.
     *  In order to do so call sec_resize(0,0);
     */ 
    void sec_resize ( const size_t size1, const size_t size2, const bool);
    
    size_t size1() const;   ///< @returns the number of rows.
    size_t size2() const;   ///< @returns the number of columns.
    size_t main_size1() const;
    size_t main_size2() const;
    size_t sec_size1() const;
    size_t sec_size2() const;
    
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
    
    void (ThisType::* row_operation_funct)( const size_t, const size_t , const size_t );
    CoefficientType& (ThisType::* op_funct)( const size_t, const size_t );
    const CoefficientType& (ThisType::* at_funct)( const size_t, const size_t ) const;
    size_t (ThisType::* size1_funct)() const;
    size_t (ThisType::* size2_funct)() const;
    
    /**
     *  In order to save we have to grad boost::serialization::access access.
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
    typedef std::vector< boost::dynamic_bitset<> > MatrixStorageType;    ///< This realizes the implementation of the data.
    typedef std::pair< size_t, size_t > MatrixEntryType;
    typedef std::list< MatrixEntryType > DiagonalType;

    MatrixBool();  ///< Creates a \f$ 0 \times 0\f$ matrix.
    /**
     *  Creates a matrix with num_rows rows and num_cols columns with entries 0.
     *  In order to get a zero matrix you may use clear();
     */
    MatrixBool( const size_t number_rows, const size_t number_cols );

    /**
     *  This performs a row operation in the Gauss algorithm.
     *  The entry in (row_1, col) is the given entry that is used to erase the entry in (row_2, col).
     */
    void row_operation( const size_t row_1, const size_t row_2, const size_t col );

    /**
     *  In order to read elements of the matrix you want to use this function.
     *  @return The function returns a copy of the given entry.
     *  @todo throw an exception if necessary i.e. if (i,j) is not a valid entry.
     */
    bool operator()( const size_t i, const size_t j );
    /**
     *  In order to keep constness, we go the usual way and implement the function at.
     *  You may want to take a look at the at()-methods of the standard containers like std::vector.
     */
    bool at( const size_t i, const size_t j ) const;

    /**
     * Adds 1 to the entry (i, j) of the matrix.
     * @note We need this function for writing access to the matrix
     * since the dynamic bitsets does not provide a reasonable reference operator.
     */
    void add_entry( const size_t i, const size_t j );

    /**
     *  As our implementation mimes ublas::matrix we use the same (awkward) method to delete a matrix.
     *  In order to do so call resize(0,0);
     */
    void resize ( const size_t size1, const size_t size2, const bool );

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

////////////////////////////////////////////////////////////////////////////////////////////////////

class MatrixBoolCSS
{
public:
    typedef bool CoefficientType;
    typedef boost::dynamic_bitset<> MatrixRowType;
    typedef std::vector< MatrixRowType > MatrixStorageType;    ///< This realizes the implementation of the data.
    typedef std::pair< size_t, size_t > MatrixEntryType;
    typedef std::list< MatrixEntryType > DiagonalType;
    typedef MatrixBoolCSS ThisType;
    
    enum MatrixBoolCSSInitialization
    {
        only_main,
        only_secondary,
        both
    };
    
    enum OperationType {
        main_and_secondary,
        secondary
    };
    
    
    MatrixBoolCSS();  ///< Creates a \f$ 0 \times 0\f$ matrix.
    
    /**
     *  Creates a matrix with num_rows rows and num_cols columns.
     *  The entries are determined by the standard constructor of CoefficientT.
     *  @warning If the standard constructor of CoefficientT does not create a coefficient with value zero (e.g. usind standard int types)
     *  the result differ from your imagination.
     *  In order to get a zero matrix you may use clear();
     */
    MatrixBoolCSS( const MatrixBoolCSSInitialization ini, const size_t num_rows1, const size_t num_cols1, const size_t num_rows2 = 0, const size_t num_cols2 = 0 );
    
    MatrixBoolCSS( const size_t number_rows, const size_t number_cols ); ///< Equivalent to MatrixBoolCSS( only_main, num_rows1, num_cols1, 0, 0 );

    /**
     *  This performs a row operation in the Gauss algorithm.
     *  The entry in (row_1, col) is the given entry that is used to erase the entry in (row_2, col).
     */
    void row_operation( const size_t row_1, const size_t row_2, const size_t col );
    void row_operation_main_and_secondary( const size_t row_1, const size_t row_2, const size_t col );
    void row_operation_secondary( const size_t row_1, const size_t row_2, const size_t col );
    void define_operations( const OperationType );

    /**
     *  In order to access elements of the matrix you want to use this function.
     *  @return The function returns a reference to the given entry.
     *  @todo throw an exception if necessary i.e. if (i,j) is not a valid entry.
     */ 
    bool operator()( const size_t i, const size_t j ); 
    bool main_op( const size_t i, const size_t j );
    bool sec_op( const size_t i, const size_t j );
    
    /**
     *  In order to keep constness, we go the usual way and implement the function at.
     *  You may want to take a look at the at()-methods of the standard containers like std::vector.
     */
    bool at( const size_t i, const size_t j ) const;
    bool main_at ( const size_t i, const size_t j) const;
    bool sec_at( const size_t i, const size_t j ) const;
    
    /**
     * Adds 1 to the entry (i, j) of the matrix.
     * @note We need this function for writing access to the matrix
     * since the dynamic bitsets does not provide a reasonable reference operator.
     */
    void add_entry( const size_t i, const size_t j );
    void main_add_entry( const size_t i, const size_t j );
    void sec_add_entry( const size_t i, const size_t j );
    
    void main_set( const size_t i, const size_t j, const bool val );
    void sec_set( const size_t i, const size_t j, const bool val );
    
    MatrixRowType& main_row( const size_t i );
    const MatrixRowType& main_row_at( const size_t i ) const;
    
    MatrixRowType& sec_row( const size_t i );
    const MatrixRowType& sec_row_at( const size_t i ) const;
    /**
     *  As our implementation mimes ublas::matrix we use the same (awkward) method to delete a matrix.
     *  In order to do so call resize(0,0);
     */ 
    void resize ( const size_t size1, const size_t size2, const bool );
    
    /**
     *  As our implementation mimes ublas::matrix we use the same (awkward) method to delete a matrix.
     *  In order to do so call sec_resize(0,0);
     */ 
    void sec_resize ( const size_t size1, const size_t size2, const bool);
    
    size_t size1() const;   ///< @returns the number of rows.
    size_t size2() const;   ///< @returns the number of columns.
    size_t main_size1() const;
    size_t main_size2() const;
    size_t sec_size1() const;
    size_t sec_size2() const;
    
    void clear();       ///< Fills every entry 0;
    void sec_clear();   ///< Fills every entry of the secondary matrix with 0;
    
     /**
      * grant std::ostream access in order to print matrix to ostreams.
     **/
    friend std::ostream& operator<< ( std::ostream& stream, const MatrixBoolCSS & matrix );

    DiagonalType diagonal;

private:
    MatrixStorageType data; ///< This realizes the data.
    MatrixStorageType sec_data; ///< This realizes the data.
    
    size_t num_rows;    ///< The number of rows.
    size_t num_cols;    ///< The number of columns.
    size_t sec_num_rows;    ///< The number of rows.
    size_t sec_num_cols;    ///< The number of columns.
    
    void (ThisType::* row_operation_funct)( const size_t, const size_t , const size_t );
    bool (ThisType::* op_funct)( const size_t, const size_t );
    bool (ThisType::* at_funct)( const size_t, const size_t ) const;
    void (ThisType::* add_entry_funct)( const size_t, const size_t );
    size_t (ThisType::* size1_funct)() const;
    size_t (ThisType::* size2_funct)() const;
    
    /**
     *  In order to save we have to grad boost::serialization::access access.
     */ 
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive &ar, const unsigned int ) ///< Implements the serialization.
    {
        ar & num_rows & num_cols & data 
           & sec_num_rows & sec_num_cols & sec_data;
    }
};

std::ostream& operator<< ( std::ostream& stream, const MatrixBoolCSS & matrix);

typedef MatrixField<Q> MatrixQ;     ///< This defines Matrices with \f$\mathbb Q\f$ coefficients.
typedef MatrixField<Zm> MatrixZm;   ///< This defines Matrices with \f$\mathbb Z/ m\mathbb Zf$ coefficients.

typedef MatrixFieldCSS<Q> MatrixCSSQ;     ///< This defines Matrices with \f$\mathbb Q\f$ coefficients.
typedef MatrixFieldCSS<Zm> MatrixCSSZm;   ///< This defines Matrices with \f$\mathbb Z/ m\mathbb Zf$ coefficients.

#include "matrix_field.ipp"

#endif // MATRIX_FIELD_HPP
