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
    typedef std::vector<CoefficientT> MatrixStorageType;    ///< This realizes the implementation of the data.
    
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
     *  As our implementation mimes ublas::matrix we use the same (awkward) method to delete a matrix.
     *  In order to do so call resize(0,0);
     */ 
    void resize (size_t size1, size_t size2, bool preserve = false)
    {
        data.clear();
        num_rows = size1;
        num_cols = size2;
        data = MatrixStorageType( size1 * size2, CoefficientT(0) );
    } 
    
    size_t size1() const;   ///< @returns the number of rows.
    size_t size2() const;   ///< @returns the number of columns.
    
    void clear();   ///< Fills every entry with CoefficientT(0);
    
     /**
     *  grant std::ostream access in order to print coefficients to ostreams like 'std::cout << Zm(44) << std::endl;'
     */
    friend std::ostream& operator<< ( std::ostream& stream, const MatrixField<CoefficientT> & matrix )
    {
        for( size_t i = 0; i < matrix.num_rows; ++i )
        {
            for( size_t j = 0; j < matrix.num_cols; )
            {
                stream << matrix.at(i,j);
                if( ++j < matrix.num_cols )
                {
                    stream << ",";
                }
                else
                {
                    stream << std::endl;
                }
            }
        }
        return stream;
    }

private:
    /**
     *  In order to keep constness, we go the usual way and implement the function at.
     *  You may want to take a look at the at()-methods of the standard containers like std::vector.
     */
    const CoefficientT& at( size_t i, size_t j ) const;
    
    MatrixStorageType data; ///< This realizes the data.
    size_t num_rows;    ///< The number of rows.
    size_t num_cols;    ///< The number of columns.
    
    /**
     *  In order to save Zm coefficients we have to grad boost::serialization::access access.
     */ 
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive &ar, const unsigned int version) ///< Implements the serialization.
    {
        ar & num_rows & num_cols & data;
    }
};

template class MatrixField<Q>;
template class MatrixField<Zm>;

typedef MatrixField<Q> MatrixQ;     ///< This defines Matrices with \f$\mathbb Q\f$ coefficients.
typedef MatrixField<Zm> MatrixZm;   ///< This defines Matrices with \f$\mathbb Z/ m\mathbb Zf$ coefficients.

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
    typedef std::vector<boost::dynamic_bitset<> > MatrixStorageType;    ///< This realizes the implementation of the data.

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
     *  In order to access elements of the matrix you want to use this function.
     *  @return The function returns a reference to the given entry.
     *  @todo throw an exception if necessary i.e. if (i,j) is not a valid entry.
     */
    bool & operator()( size_t i, size_t j );

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
    void resize (size_t size1, size_t size2, bool preserve = false)
    {
        data.clear();
        std::vector< boost::dynamic_bitset<> > new_data( size1, boost::dynamic_bitset<>(size2, 0) );
        std::swap(data, new_data);
    }

    size_t size1() const;   ///< @returns the number of rows.
    size_t size2() const;   ///< @returns the number of columns.

    void clear();   ///< Fills every entry with 0;

     /**
     *  grant std::ostream access in order to print coefficients to ostreams like 'std::cout << Zm(44) << std::endl;'
     */
    friend std::ostream& operator<< ( std::ostream& stream, const MatrixBool & matrix)
    {
        for( size_t i = 0; i < matrix.num_rows; ++i )
        {
            for( size_t j = 0; j < matrix.num_cols; )
            {
                stream << matrix.at(i,j);
                if( ++j < matrix.num_cols )
                {
                    stream << ",";
                }
                else
                {
                    stream << std::endl;
                }
            }
        }
        return stream;
    }

private:
    /**
     *  In order to keep constness, we go the usual way and implement the function at.
     *  You may want to take a look at the at()-methods of the standard containers like std::vector.
     */
    const bool & at( size_t i, size_t j ) const;

    MatrixStorageType data; ///< This realizes the data.
    size_t num_rows;    ///< The number of rows.
    size_t num_cols;    ///< The number of columns.

    /**
     *  In order to save coefficients we have to grad boost::serialization::access access.
     */
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive &ar, const unsigned int version) ///< Implements the serialization.
    {
        ar & num_rows & num_cols & data;
    }
};

#include "matrix_field.ipp"
#endif // MATRIX_FIELD_HPP
