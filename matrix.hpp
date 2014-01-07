#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <vector>

#include "parallelization.hpp"

template < class CoefficientT >
class Matrix
{
public:
    typedef std::vector<CoefficientT> MatrixStorageType;
    
    Matrix();
    Matrix( size_t number_rows, size_t number_cols );
    
    void row_operation( size_t row_1, size_t row_2, size_t col );
    
    CoefficientT & operator()( size_t i, size_t j );

    size_t size1() const;
    size_t size2() const;
    
    friend std::ostream& operator<< ( std::ostream& stream, const Matrix<CoefficientT> & matrix )
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
    const CoefficientT& at( size_t i, size_t j ) const;
    
    MatrixStorageType data;
    size_t num_rows;
    size_t num_cols;
    
    friend class boost::serialization::access;
    
    template <class Archive>
    void serialize(Archive &ar, const unsigned int version) ///< Implements the serialization.
    {
        ar & num_rows & num_cols & data;
    }
};

#include "matrix.ipp"
#endif // MATRIX_HPP
