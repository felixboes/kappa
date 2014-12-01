#include "diagonalizer_field.hpp"

template< class MatrixType >
uint32_t DiagonalizerField< MatrixType >::dfct() const
{
    return def;
}

template< class MatrixType >
HomologyField::KernT DiagonalizerField< MatrixType >::kern()
{
    return def;
}

template< class MatrixType >
uint32_t DiagonalizerField< MatrixType >::rank() const
{
    return rnk;
}

template< class MatrixType >
HomologyField::TorsT DiagonalizerField< MatrixType >::tors()
{
    return rnk;
}

template< class MatrixType >
void DiagonalizerField< MatrixType > :: apply_base_changes( MatrixType& differential, const MatrixType& base_changes )
{
    if( base_changes.size1() == 0 ||
        base_changes.size2() == 0 ||
        differential.size1() == 0 ||
        differential.size2() == 0
      )
    {
        return;
    }
    
    
    if( transp == false )
    {
        const auto num_rows = differential.size1();
        const auto num_cols = differential.size2();
        if( base_changes.size2() != num_rows )
        {
            std::cout << "Error: The base changes have " << base_changes.size2() << " many columns but the current differential has " << differential.size1() << " many rows." << std::endl;
            return;
        }
        
        // The column operations can be read off the triangular shape.
        VectorField< CoefficientType > base_change_single_col( num_rows );
        
        for( auto diag_entry : base_changes.diagonal )
        {
            // all entries left of diag_entry are zero.
            // The k-th row of the differential is altered using the entries right of the diagonal entry.
            const size_t k = diag_entry.first;
            const size_t l = diag_entry.second;
            const CoefficientType lambda = CoefficientType(1) / base_changes.at( k, l );
            std::vector< Thread > threads;
            
            // Prepare one of the base changes.
            // Parallelized version of:
            //   for( size_t i = l+1; i < num_rows; ++i )
            //   {
            //       base_change_single_col(i) = lambda * base_changes.at( k, i );
            //   }
            size_t chunk_size = ( num_rows - (l+1) ) / num_working_threads;
            if( ( num_rows - (l+1) ) % num_working_threads != 0 )
            {
                chunk_size++;
            }
            
            bool last_round = false;
            for( size_t i = 0; last_round == false; ++i )
            {
                const size_t offset = l+1 + i*chunk_size;
                // recompute chunk_size if nessecary.
                if( chunk_size*(i+1) >= num_rows - (l+1) )
                {
                    chunk_size = num_rows - (l+1) - i*chunk_size;
                    last_round = true;
                }
                
                threads.emplace_back(
                    [&base_change_single_col, &base_changes, k, lambda, offset, chunk_size]
                    {
                        for( size_t j = 0; j < chunk_size; ++j )
                        {
                            base_change_single_col(offset + j) = lambda * base_changes.at( k, offset + j );
                        }
                    }
                );
            }
            
            for( auto& it : threads )
            {
                it.execute();
            }
            for( auto& it : threads )
            {
                it.wait();
            }
            threads.clear();
            
            // Apply base change.
            // Parallized version of:
            //
            // for( size_t j = 0; j < num_cols; ++j )
            // {
            //     CoefficientType& altered_entry = differential( l, j );
            //     for( size_t i = l+1; i < num_rows; ++i )
            //     {
            //         altered_entry += base_change_single_col(i) * differential.at(i,j);
            //     }
            //  }
            chunk_size = num_cols / num_working_threads;
            if( num_cols  % num_working_threads != 0 )
            {
                chunk_size++;
            }
            
            last_round = false;
            for( size_t i = 0; last_round == false; ++i )
            {
                const size_t offset = i*chunk_size;
                // recompute chunk_size if nessecary.
                if( chunk_size*(i+1) >= num_cols )
                {
                    chunk_size = num_cols - i*chunk_size;
                    last_round = true;
                }
                
                threads.emplace_back(
                    [&differential, &base_change_single_col, l, lambda, offset, chunk_size, num_rows]
                    {
                        for( size_t j = 0; j < chunk_size; ++j )
                        {
                            for( size_t i = l+1; i < num_rows; ++i )
                            {
                                differential( l, offset + j ) += base_change_single_col( i ) * differential.at( i, offset + j );
                            }
                        }
                    }
                );
            }
            
            for( auto& it : threads )
            {
                it.execute();
            }
            for( auto& it : threads )
            {
                it.wait();
            }
            threads.clear();
        }
    }
    else
    {
        std::cout << "Todo." << std::endl;
        
    }
}


/**
 *  Compute the dimension of the image of matrix.
 *  This done by computing the number of lineary independent columns or rows.
 */
template< class MatrixType >
uint32_t DiagonalizerField< MatrixType >::diag_field( MatrixType &matrix )
{
    size_t num_rows = matrix.size1();
    size_t num_cols = matrix.size2();

    /// @TODO: Is #rows <= #cols?
    // if( num_rows <= num_cols )
    {
        // We avoid permutations of rows in the Gauss algorithm by performing the following steps:
        //   - Iterate through the columns.
        //   - Iterate through a singly linked vector of rows.
        //   - For a given column, remove a row from the vector if the common entry of column and row is invertible.
        //   - Apply row operations to all rows beneath the row.

        std::vector<size_t> rows_to_check;

        // Fill the vector.
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
                matrix.diagonal.emplace_back( row_1, col );
                // We do not need to check this row in further iterations.
                // After erasing, it will point to the next element in the vector.
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
void DiagonalizerField< MatrixType >::operator() ( MatrixType& matrix )
{
    current_rank = 0;
    if( num_working_threads == 1 && num_remaining_threads == 0 )
    {
        rnk = diag_field( matrix );
    }
    else
    {
        // Parallelized diagonalization needs at least one remaining_thread.
        if (num_remaining_threads == 0)
        {
            --num_working_threads;
            ++num_remaining_threads;
        }
        rnk = diag_field_parallelized( matrix );
    }
    if( transp == false )
    {
        def = matrix.size2() - ommit_rows.size() - rnk;
    }
    else
    {
        def = matrix.size1() - ommit_rows.size() - rnk;
    }
}


////////////////////////////////////////////////////////////////////////////////


template < class MatrixType >
DiagonalizerField< MatrixType >::JobQueue::JobQueue(
    MatrixType&                         matrix_init,
    typename MatrixType::DiagonalType&  diag,
    const std::list< size_t >&          ommit_rows,
    const uint32_t                      number_of_working_threads,
    const uint32_t                      number_of_remaining_threads,
    atomic_uint&                        current_rank )
:
    matrix(matrix_init),
    rows_to_work_at(),
    remaining_rows(),
    work_chunk_size(0),
    remaining_chunk_size(0),
    num_working_threads(number_of_working_threads),
    num_remaining_threads(number_of_remaining_threads),
    row_1(0),
    col(0),
    diagonal(diag)
{
    if ( matrix.size1() == 0 || matrix.size2() == 0 ) // Then no diagonalizing is necessary.
    {
        return;
    }
    // Sort the entries of the first column by whether they are zero or not.
    auto ommit_it = ommit_rows.begin();
    for( size_t row = 0 ; row < matrix.size1(); ++row )
    {
        if( ommit_it != ommit_rows.end() && row == *ommit_it )
        {
            ++ommit_it;
            continue;
        }
        if ( matrix(row, col) == typename MatrixType::CoefficientType(0) )
        {
            remaining_rows.push_back(row);
        }
        else
        {
            rows_to_work_at.push_back(row);
        }
    }

    // If we have found a row with non-zero entry in the first column, take
    // one row as the row which is added to all the other rows in the set of
    // row operations for the first column.
    if (rows_to_work_at.size() != 0)
    {
        row_1 = rows_to_work_at.back();
        diagonal.emplace_back( row_1, col );
        rows_to_work_at.pop_back();
        // This row contributes to the rank.
        ++current_rank;
    }

    recompute_chunk_sizes();
}

#ifdef BROKEN_VECTOR_IMPLEMENTATION
template < class MatrixType >
DiagonalizerField< MatrixType >::JobQueue::JobQueue(JobQueue const & other)
:
   matrix(other.matrix),
   rows_to_work_at(std::move(other.rows_to_work_at)),
   remaining_rows(std::move(other.remaining_rows)),
   work_chunk_size(other.work_chunk_size),
   remaining_chunk_size(other.remaining_chunk_size),
   num_working_threads(other.num_working_threads),
   num_remaining_threads(other.num_remaining_threads),
   row_1(other.row_1),
   col(other.col)
{}

template < class MatrixType >
void DiagonalizerField< MatrixType >::JobQueue::operator=(JobQueue const & other)
{
   matrix                 = other.matrix;
   rows_to_work_at        = other.rows_to_work_at;
   remaining_rows         = other.remaining_rows;
   work_chunk_size        = other.work_chunk_size;
   remaining_chunk_size   = other.remaining_chunk_size;
   num_working_threads    = other.num_working_threads;
   num_remaining_threads  = other.num_remaining_threads;
   row_1                  = other.row_1;
   col                    = other.col;
}
#endif

/**
 * Clears the vectors new_rows_to_work_at and new_remaining_rows from each worker
 * and returns their union in one vector each.
 */
template < class MatrixType >
void DiagonalizerField< MatrixType >::JobQueue:: add_up_rows
    (std::vector<Worker> & workers,
     std::vector<size_t> & former_rows_to_work_at,
     std::vector<size_t> & former_remaining_rows)
{
    std::vector<size_t> rows_to_work_at;
    std::vector<size_t> remaining_rows;

    // Determine the total size of the vectors of each worker.
    size_t total_num_rows_to_work_at = 0;
    size_t total_num_remaining_rows  = 0;

    for (auto & it : workers)
    {
        Worker & worker = it;
        total_num_rows_to_work_at += worker.get_new_rows_to_work_at().size();
        total_num_remaining_rows  += worker.get_new_remaining_rows().size();
    }

    // Reserve memory for the union of all vectors.
    rows_to_work_at.reserve(total_num_rows_to_work_at);
    remaining_rows.reserve (total_num_remaining_rows);

    for (auto & it : workers)
    {
        Worker & worker = it ;
        std::vector<size_t> & new_rows_to_work_at = worker.get_new_rows_to_work_at();
        std::vector<size_t> & new_remaining_rows  = worker.get_new_remaining_rows();
        // Append the vectors of this worker to the united array.
        rows_to_work_at.insert( rows_to_work_at.end(), new_rows_to_work_at.begin(), new_rows_to_work_at.end() );
        remaining_rows.insert ( remaining_rows.end(),  new_remaining_rows.begin(),  new_remaining_rows.end() );
        // Clear the vectors of this worker.
        new_rows_to_work_at.clear();
        new_remaining_rows.clear();
    }

    rows_to_work_at.swap(former_rows_to_work_at);
    remaining_rows.swap(former_remaining_rows);
}

template < class MatrixType >
bool DiagonalizerField< MatrixType >::JobQueue::update_rank_and_work( atomic_uint & current_rank, std::vector< Worker > & workers )
{
    ++col;

    if ( col < matrix.size2() && current_rank < matrix.size1() )
    {
        add_up_rows( workers, rows_to_work_at, remaining_rows );
        if ( rows_to_work_at.size() != 0 )
        {
            row_1 = rows_to_work_at.back();
            diagonal.push_back( MatrixEntryType(row_1, col) );
            rows_to_work_at.pop_back();
            ++current_rank;
        }

        recompute_chunk_sizes();

        return true;
    }
    else
    {
        return false;
    }

}

template < class MatrixType >
void DiagonalizerField< MatrixType >::JobQueue::get_work_chunk(
    uint32_t const  thread_id,
    size_t &        begin,
    size_t &        end) const
{
    begin = std::min(  thread_id      * work_chunk_size, rows_to_work_at.size() );
    end   = std::min( (thread_id + 1) * work_chunk_size, rows_to_work_at.size() );
}

template < class MatrixType >
void DiagonalizerField< MatrixType >::JobQueue::get_remaining_chunk(
    uint32_t const  thread_id,
    size_t &        begin,
    size_t &        end) const
{
    begin = std::min(  thread_id      * remaining_chunk_size, remaining_rows.size() );
    end   = std::min( (thread_id + 1) * remaining_chunk_size, remaining_rows.size() );
}

template < class MatrixType >
void DiagonalizerField< MatrixType >::JobQueue::recompute_chunk_sizes()
{
    work_chunk_size = rows_to_work_at.size() / (num_working_threads);
    if ( rows_to_work_at.size() % (num_working_threads) != 0 )
    {
        ++work_chunk_size;
    }

    remaining_chunk_size = remaining_rows.size() / (num_remaining_threads);
    if ( remaining_rows.size() % num_remaining_threads != 0 )
    {
        ++remaining_chunk_size;
    }
}


////////////////////////////////////////////////////////////////////////////////


template< class MatrixType >
DiagonalizerField< MatrixType >::Worker::Worker( const uint32_t identification, JobQueue &sl )
:
    id(identification),
    jobs(sl),
    new_rows_to_work_at(),
    new_remaining_rows()
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
void DiagonalizerField< MatrixType >::Worker::work( MatrixType& matrix )
{
    size_t begin = 0, end = 0;
    size_t new_col = jobs.col + 1;
    for ( jobs.get_work_chunk(id, begin, end); begin < end; ++begin )
    {
        size_t row_2 = jobs.rows_to_work_at[begin];
        matrix.row_operation( jobs.row_1, row_2, jobs.col );
    }
    if (new_col < matrix.size2())
    {
        for ( jobs.get_work_chunk(id, begin, end); begin < end; ++begin )
        {
            size_t row_2 = jobs.rows_to_work_at[begin];
            if ( matrix(row_2, new_col) == typename MatrixType::CoefficientType(0) )
            {
                new_remaining_rows.push_back(row_2);
            }
            else
            {
                new_rows_to_work_at.push_back(row_2);
            }
        }
    }
}

template< class MatrixType >
std::vector<size_t> & DiagonalizerField< MatrixType >::Worker::get_new_rows_to_work_at()
{
    return new_rows_to_work_at;
}

template< class MatrixType >
std::vector<size_t> & DiagonalizerField< MatrixType >::Worker::get_new_remaining_rows()
{
    return new_remaining_rows;
}

template< class MatrixType >
void DiagonalizerField< MatrixType >::Worker::collect_remaining_work( MatrixType & matrix )
{
    std::vector<size_t> & remaining_rows = jobs.remaining_rows;
    size_t begin = 0;
    size_t end = 0;
    size_t new_col = jobs.col + 1;
    if ( new_col < matrix.size2() )
    {
        for ( jobs.get_remaining_chunk(id, begin, end); begin < end; ++begin )
        {
            size_t row = remaining_rows[begin];
            if ( matrix(row, new_col) == typename MatrixType::CoefficientType(0) )
            {
                new_remaining_rows.push_back(row);
            }
            else
            {
                new_rows_to_work_at.push_back(row);
            }
        }
    }
}

template< class MatrixType >
uint32_t DiagonalizerField< MatrixType >::diag_field_parallelized( MatrixType & matrix )
{
    JobQueue jobs(matrix, matrix.diagonal, ommit_rows, num_working_threads, num_remaining_threads, current_rank);
    std::vector<Worker> workers;
    std::vector<Thread> threads;

    // Initialize num_working_threads threads to perform row operations.
    for( size_t i = 0; i < num_working_threads; ++i )
    {
        workers.emplace_back( i, jobs );
        threads.emplace_back( [i, &workers, &matrix]{workers[i].work(matrix);} );
    }

    // Initialize num_remaining_threads further threads to check for row operations of the next column
    // among the rows for which no row operation is performed.
    for ( size_t i = 0; i < num_remaining_threads; ++i)
    {
        workers.emplace_back(i, jobs);
        size_t thread_id = num_working_threads + i;
        threads.emplace_back( [thread_id, &workers, &matrix]{workers[thread_id].collect_remaining_work(matrix);} );
    }
    // Start threads.
    do
    {
        for ( auto & it : threads )
        {
            it.execute();
        }
        // Wait to finish.
        for( auto & it : threads )
        {
            it.wait();
        };
    } while ( jobs.update_rank_and_work( current_rank, std::ref(workers) ) );
    return current_rank;
}
