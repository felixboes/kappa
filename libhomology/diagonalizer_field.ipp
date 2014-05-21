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
void DiagonalizerField< MatrixType >::operator() ( MatrixType &matrix, uint32_t number_threads )
{
    // call operator() ( MatrixType &m, atomic_uint &, uint32_t ).
    atomic_uint current_rank;
    this->operator ()(matrix, current_rank, number_threads);
}

template< class MatrixType >
void DiagonalizerField< MatrixType >::operator() ( MatrixType &matrix, atomic_uint & current_rank, uint32_t num_working_threads, uint32_t num_remaining_threads )
{
    // We use the non-parallel version.
    if( num_working_threads == 1 && num_remaining_threads == 0)
    {
        rnk = diag_field(matrix, current_rank);
    }
    else
    {
        // Parallelized diagonalization needs at least one remaining_thread.
        if (num_remaining_threads == 0)
        {
            --num_working_threads;
            ++num_remaining_threads;
        }
        rnk = diag_field_parallelized( matrix, current_rank, num_working_threads, num_remaining_threads );
    }
    def = matrix.size1() - rnk;
}


////////////////////////////////////////////////////////////////////////////////


template < class MatrixType >
DiagonalizerField< MatrixType >::JobQueue::JobQueue(MatrixType & matrix_init, uint32_t number_of_working_threads, uint32_t number_of_remaining_threads, atomic_uint & current_rank )
:
    matrix(matrix_init),
    rows_to_work_at(),
    remaining_rows(),
    work_chunk_size(0),
    remaining_chunk_size(0),
    num_working_threads(number_of_working_threads),
    num_remaining_threads(number_of_remaining_threads),
    row_1(0),
    col(0)
{
    if (matrix.size2() == 0) // Then no diagonalizing is necessary.
    {
        return;
    }
    // Sort the entries of the first column by whether they are zero or not.
    for( size_t row = 0; row < matrix.size1(); ++row )
    {
        if ( matrix(row, col) == typename MatrixType::CoefficientType(0))
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
        rows_to_work_at.insert(rows_to_work_at.end(), new_rows_to_work_at.begin(), new_rows_to_work_at.end());
        remaining_rows.insert (remaining_rows.end(),  new_remaining_rows.begin(),  new_remaining_rows.end());
        // Clear the vectors of this worker.
        new_rows_to_work_at.clear();
        new_remaining_rows.clear();
    }

    rows_to_work_at.swap(former_rows_to_work_at);
    remaining_rows.swap(former_remaining_rows);
}

template < class MatrixType >
bool DiagonalizerField< MatrixType >::JobQueue::update_rank_and_work(atomic_uint & current_rank, std::vector<Worker> & workers)
{
    ++col;

    if (col < matrix.size2() && current_rank < matrix.size1())
    {
        add_up_rows(workers, rows_to_work_at, remaining_rows);
        if (rows_to_work_at.size() != 0)
        {
            row_1 = rows_to_work_at.back();
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
    begin = std::min( thread_id      * work_chunk_size, rows_to_work_at.size());
    end   = std::min((thread_id + 1) * work_chunk_size, rows_to_work_at.size());
}

template < class MatrixType >
void DiagonalizerField< MatrixType >::JobQueue::get_remaining_chunk(
    uint32_t const  thread_id,
    size_t &        begin,
    size_t &        end) const
{
    begin = std::min( thread_id      * remaining_chunk_size, remaining_rows.size());
    end   = std::min((thread_id + 1) * remaining_chunk_size, remaining_rows.size());
}

template < class MatrixType >
void DiagonalizerField< MatrixType >::JobQueue::recompute_chunk_sizes()
{
    work_chunk_size = rows_to_work_at.size() / (num_working_threads);
    if (rows_to_work_at.size() % (num_working_threads) != 0)
    {
        ++work_chunk_size;
    }

    remaining_chunk_size = remaining_rows.size() / (num_remaining_threads);
    if (remaining_rows.size() % num_remaining_threads != 0)
    {
        ++remaining_chunk_size;
    }
}


////////////////////////////////////////////////////////////////////////////////


template< class MatrixType >
DiagonalizerField< MatrixType >::Worker::Worker(uint32_t identification, JobQueue &sl)
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
void DiagonalizerField< MatrixType >::Worker::work(MatrixType& matrix)
{
    size_t begin = 0, end = 0;
    size_t new_col = jobs.col + 1;
    bool compute_next_row = (new_col < matrix.size2());
    for ( jobs.get_work_chunk(id, begin, end); begin < end; ++begin )
    {
        size_t row_2 = jobs.rows_to_work_at[begin];
        matrix.row_operation( jobs.row_1, row_2, jobs.col );
        // ToDo: Improve performance!
        if (compute_next_row)
        {
            if (matrix(row_2, new_col) == typename MatrixType::CoefficientType(0))
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
void DiagonalizerField< MatrixType >::Worker::collect_remaining_work(MatrixType & matrix)
{
    std::vector<size_t> & remaining_rows = jobs.remaining_rows;
    size_t begin = 0;
    size_t end = 0;
    size_t new_col = jobs.col + 1;
    if (new_col < matrix.size2())
    {
        for (jobs.get_remaining_chunk(id, begin, end); begin < end; ++begin)
        {
            size_t row = remaining_rows[begin];
            if (matrix(row, new_col) == typename MatrixType::CoefficientType(0))
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
uint32_t DiagonalizerField< MatrixType >::diag_field_parallelized(
    MatrixType &   matrix,
    atomic_uint &  current_rank,
    uint32_t       num_working_threads,
    uint32_t       num_remaining_threads)
{
    JobQueue jobs(matrix, num_working_threads, num_remaining_threads, current_rank);
    std::vector<Worker> workers;
    std::vector<Thread> threads;

    // Initialize num_working_threads threads to perform row operations.
    for( size_t i = 0; i < num_working_threads; ++i )
    {
        workers.emplace_back( i, jobs );
        threads.emplace_back( [i, &matrix, &workers]{workers[i].work(matrix);} );
    }

    // Initialize num_remaining_threads further threads to check for row operations of the next column
    // among the rows for which no row operation is performed.
    for ( size_t i = 0; i < num_remaining_threads; ++i)
    {
        workers.emplace_back(i, jobs);
        threads.emplace_back([i, num_working_threads, &matrix, &workers]{workers[num_working_threads + i].collect_remaining_work(matrix);});
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
    }while (jobs.update_rank_and_work(current_rank, std::ref(workers)));
    return current_rank;
}
