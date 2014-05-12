#include "diagonalizer_field.hpp"

template< class MatrixType >
uint32_t DiagonalizerField< MatrixType >::dfct()
{
    return def;
}

template< class MatrixType >
HomologyField::KernT DiagonalizerField< MatrixType >::kern()
{
    return def;
}

template< class MatrixType >
uint32_t DiagonalizerField< MatrixType >::rank()
{
    return rnk;
}

template< class MatrixType >
HomologyField::TorsT DiagonalizerField< MatrixType >::tors()
{
    return rnk;
}

/**
 *  Compute the dimension of the image of matrix.
 *  This done by computing the number of lineary independent columns or rows.
 */
template< class MatrixType >
uint32_t DiagonalizerField< MatrixType >::diag_field(MatrixType &matrix, atomic_uint& current_rank)
{
    size_t num_rows = matrix.size1();
    size_t num_cols = matrix.size2();
    current_rank = 0; // Stores the dimension of the subspace spanned by the first 0 \le j < col columns. The dimension is at most of size row.

    /// @TODO: Is #rows <= #cols?
    // if( num_rows <= num_cols )
    {
        // We avoid permutations of rows in the Gauss algorithm by performing the following steps:
        //   - Iterate through the columns.
        //   - Iterate through a singly linked list of rows.
        //   - For a given column, remove a row from the list if the common entry of column and row is invertible.
        //   - Apply row operations to all rows beneath the row.

        std::list<size_t> rows_to_check;

        // Fill the list.
        for( size_t row_1 = 0; row_1 < num_rows; ++row_1 )
        {
            rows_to_check.push_back(row_1);
        }

        size_t col = 0;
        size_t row_1 = 0;

        // Iterate through columns. We may stop if the rank is maximal.
        for (; col < num_cols && current_rank < num_rows; ++col )
        {
            // Find first invertible (i.e. non-zero) entry in the remaining rows.
            auto it = rows_to_check.begin();
            for( ; it != rows_to_check.end() && matrix( *it, col ) == typename MatrixType::CoefficientType(0); ++it )
            {
            }

            // Does this row contribute to the rank, i.e. did we find an invertible element?
            if( it != rows_to_check.end() )
            {
                // This row contributes to the rank.
                current_rank++;
                row_1 = *it;
                // We do not need to check this row in further iterations.
                // After erasing, it will point to the next element in the list.
                it = rows_to_check.erase(it);
                // Use row operations to zeroize the remaining elements of the column.
                for( size_t relative_position=0 ; it != rows_to_check.end(); ++it, ++relative_position )
                {
                    size_t row_2 = *it;
                    // Assuming that the entry (row_2, col) differs from zero, we perform
                    // a row operation to matrix to zeroize the entry (row_2, col) using the entry (row_1, col).
                    if( matrix( row_2, col ) != typename MatrixType::CoefficientType(0) )
                    {
                        matrix.row_operation( row_1, row_2, col );
                    }
                }
            }
        }
    }
    return current_rank;
}

template< class MatrixType >
void DiagonalizerField< MatrixType >::operator() ( MatrixType &matrix, uint32_t number_threads )
{
    // call operator() ( MatrixType &m, atomic_uint &, uint32_t ).
    atomic_uint current_rank;
    this->operator ()(matrix, current_rank, number_threads);
}

template< class MatrixType >
void DiagonalizerField< MatrixType >::operator() ( MatrixType &matrix, atomic_uint & current_rank, uint32_t number_threads )
{
    if( number_threads == 0 )
    {
        rnk = diag_field(matrix, current_rank);
    }
    else
    {
        rnk = diag_field_parallelized( matrix, current_rank, number_threads );
    }
    def = matrix.size1() - rnk;
}


////////////////////////////////////////////////////////////////////////////////


template < class MatrixType >
DiagonalizerField< MatrixType >::JobQueue::JobQueue(MatrixType & matrix_init, uint32_t number_of_working_threads )
:
    matrix(matrix_init),
    number_of_threads(number_of_working_threads),
    rows_to_check(matrix.size1()),
    col(0)
{
    for( size_t row = 0; row < matrix.size1(); ++row )
    {
        rows_to_check[row] = row;
    }
}

#ifdef BROKEN_VECTOR_IMPLEMENTATION
template < class MatrixType >
DiagonalizerField< MatrixType >::JobQueue::JobQueue(JobQueue const & other)
:
   matrix(other.matrix),
   number_of_threads(other.number_of_threads),
   rows_to_check(std::move(other.rows_to_check)),
   col(other.col)
{}

template < class MatrixType >
void DiagonalizerField< MatrixType >::JobQueue::operator=(JobQueue const & other)
{
   matrix             = other.matrix;
   number_of_threads  = other.number_of_threads;
   rows_to_check      = other.rows_to_check;
   col                = other.col;
}
#endif



template < class MatrixType >
typename DiagonalizerField< MatrixType >::RowOpParam const &
DiagonalizerField< MatrixType >::JobQueue::get_operation(size_t op_id) const
{
    return operations[op_id];
}

template < class MatrixType >
void DiagonalizerField< MatrixType >::JobQueue::get_chunk(
    uint32_t const  thread_id,
    size_t &        begin,
    size_t &        end) const
{
    begin = std::min( thread_id      * cur_chunk_size, operations.size());
    end   = std::min((thread_id + 1) * cur_chunk_size, operations.size());
}

template< class MatrixType >
bool DiagonalizerField< MatrixType >::JobQueue::compute_operations(atomic_uint & current_rank)
{
    operations.clear();

    size_t num_rows = matrix.size1();
    size_t num_cols = matrix.size2();

    // Iterate through columns. We may stop if the rank is maximal.
    for (; col < num_cols && current_rank < num_rows; ++col )
    {

        // Find first invertible (i.e. non-zero) entry in the remaining rows.
        auto it = rows_to_check.begin();
        for( ; it != rows_to_check.end() && matrix( *it, col ) == typename MatrixType::CoefficientType(0); ++it )
        {
        }

        // Does this row contribute to the rank, i.e. did we find an invertible element?
        if( it != rows_to_check.end() )
        {
            // This row contributes to the rank.
            current_rank++;
            size_t row_1 = *it;
            // We do not need to check this row in further iterations.
            // After erasing, it will point to the next element in the list.
            it = rows_to_check.erase(it);
            // Use row operations to zeroize the remaining elements of the column.
            for( size_t relative_position=0 ; it != rows_to_check.end(); ++it, ++relative_position )
            {
                size_t row_2 = *it;
                // Assuming that the entry (row_2, col) differs from zero, we perform
                // a row operation to matrix to zeroize the entry (row_2, col) using the entry (row_1, col).
                if( matrix( row_2, col ) != typename MatrixType::CoefficientType(0) )
                {
                    operations.emplace_back( row_1, row_2, col );
                }
            }
            ++col;
            recompute_chunk_size();
            return true;
        }
    }
    return false;
}

template < class MatrixType >
void DiagonalizerField< MatrixType >::JobQueue::recompute_chunk_size()
{
    cur_chunk_size = operations.size() / number_of_threads;
    if (operations.size() % number_of_threads != 0)
    {
        ++cur_chunk_size;
    }
}


////////////////////////////////////////////////////////////////////////////////


template< class MatrixType >
DiagonalizerField< MatrixType >::Worker::Worker(uint32_t identification, JobQueue &sl)
:
    id(identification),
    jobs(sl)
{}

#ifdef BROKEN_VECTOR_IMPLEMENTATION
template< class MatrixType >
DiagonalizerField< MatrixType >::Worker::Worker(Worker const & other)
:
    id(other.id),
    jobs(other.jobs)
{}

template< class MatrixType >
void DiagonalizerField< MatrixType >::Worker::operator=(Worker const & other)
{
    id = other.id;
    jobs = other.jobs;
}
#endif

template< class MatrixType >
void DiagonalizerField< MatrixType >::Worker::work(MatrixType& matrix)
{
    typedef DiagonalizerField<MatrixType>::RowOpParam Operation;

    size_t begin = 0, end = 0;
    for ( jobs.get_chunk(id, begin, end); begin < end; ++begin )
    {
        Operation const & op = jobs.get_operation(begin);
        matrix.row_operation( op.row1, op.row2, op.col );
    }
}

template< class MatrixType >
uint32_t DiagonalizerField< MatrixType >::diag_field_parallelized(
    MatrixType &   matrix,
    atomic_uint &  current_rank,
    uint32_t       number_threads )
{
    JobQueue jobs(matrix, number_threads);
    std::vector<Worker> workers;
    std::vector<Thread> threads;
    for( size_t i = 0; i < number_threads; ++i )
    {
        workers.emplace_back( i, jobs );
        threads.emplace_back( [i, &matrix, &workers]{workers[i].work(matrix);} );
    }

    // Start threads.
    current_rank = 0;
    while (jobs.compute_operations(current_rank))
    {
        for ( auto & it : threads )
        {
            it.execute();
        }

        // Wait to finish.
        for( auto & it : threads )
        {
            it.wait();
        }
    }
    return current_rank;
}
