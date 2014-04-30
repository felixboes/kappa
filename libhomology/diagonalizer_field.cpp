#include "diagonalizer_field.hpp"

uint32_t DiagonalizerBool::dfct()
{
    return def;
}

HomologyField::KernT DiagonalizerBool::kern()
{
    return def;
}

uint32_t DiagonalizerBool::rank()
{
    return rnk;
}

HomologyField::TorsT DiagonalizerBool::tors()
{
    return rnk;
}

/**
 *  Compute the dimension of the image of matrix.
 *  This done by computing the number of lineary independent columns or rows.
 */
uint32_t DiagonalizerBool::diag_field(MatrixBool &matrix, atomic_uint& current_rank)
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
            for( ; it != rows_to_check.end() && matrix( *it, col ) == 0; ++it )
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
                    if( matrix( row_2, col ) != 0 )
                    {
                        matrix.row_operation( row_1, row_2, col );
                    }
                }
            }
        }
    }
    return current_rank;
}

void DiagonalizerBool::operator() ( MatrixBool &matrix, uint32_t number_threads )
{
    // call operator() ( MatrixBool &m, atomic_uint &, uint32_t ).
    atomic_uint current_rank;
    this->operator ()(matrix, current_rank, number_threads);
}

void DiagonalizerBool::operator() ( MatrixBool &matrix, atomic_uint & current_rank, uint32_t number_threads )
{
    if( number_threads == 0 )
    {
        rnk = diag_field(matrix, current_rank);
    }
    else
    {
        rnk = diag_field_parallelized( matrix, current_rank, number_threads );
    }
    def = matrix.size2() - rnk;
}

void DiagonalizerBool::SyncList::all_work_done()
{
    std::lock_guard<std::mutex> lock(mtx);
    work_done = true;
    cond_wait_for_filler.notify_all();
}

bool DiagonalizerBool::SyncList::no_work_left()
{
    std::lock_guard<std::mutex> lock(mtx);
    return work_done == true;
}

void DiagonalizerBool::SyncList::put( DiagonalizerBool::RowOpParam rop )
{
    std::lock_guard<std::mutex> lock(mtx);
    row_op_list.push_back(rop);
    cond_wait_for_filler.notify_one();
}

bool DiagonalizerBool::SyncList::get( DiagonalizerBool::RowOpParam & rop )
{
    std::unique_lock<std::mutex> lock(mtx);
    // Notify filler thread if necessary.
    if( row_op_list.empty() == true )
    {
        // Worker is not busy anymore.
        num_busy_workers--;
        cond_wait_for_worker.notify_one();

        // wait blocks as long the lambda function returns false, i.e. as long as the condition A in {return !( A );} is true.
        cond_wait_for_filler.wait( lock, [&]{ return ! (row_op_list.empty() == true && work_done == false );} );

        // Worker will obtain the next job.
        num_busy_workers++;
    }

    // Why did we continue? Is there still work?
    if( work_done == false )
    {
        rop = row_op_list.front();
        row_op_list.pop_front();
        return true;
    }
    else
    {
        return false;
    }
}

void DiagonalizerBool::SyncList::wait_for_workers()
{
    std::unique_lock<std::mutex> lock(mtx);
    // Wait till the list is empty.
    cond_wait_for_worker.wait( lock, [&]{ return !( row_op_list.empty() == false ); } );

    // Wait for worker threads to finish their jobs.
    cond_wait_for_worker.wait( lock, [&]{ return !( num_busy_workers > 0 ); } );
}

DiagonalizerBool::Filler::Filler( SyncList &sl ) : sync_list(sl) {}

void DiagonalizerBool::Filler::work(MatrixBool& matrix, atomic_uint& current_rank)
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
            for( ; it != rows_to_check.end() && matrix( *it, col ) == 0; ++it )
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
                    if( matrix( row_2, col ) != 0 )
                    {
                        sync_list.put( DiagonalizerBool::RowOpParam( row_1, row_2, col ) );
                    }
                }
                // Wait till all work for the current column is done.
                sync_list.wait_for_workers();
            }
        }
    }

    // Tell worker threads that there will be no more work.
    sync_list.all_work_done();
}

DiagonalizerBool::Worker::Worker(uint32_t identification, SyncList &sl) : id(identification), sync_list(sl) {}

void DiagonalizerBool::Worker::work(MatrixBool& matrix)
{
    DiagonalizerBool::RowOpParam rop;

    while( sync_list.get(rop) == true )
    {
        matrix.row_operation(rop.row1, rop.row2, rop.col);
    }
}

uint32_t DiagonalizerBool::diag_field_parallelized(MatrixBool& matrix, atomic_uint & current_rank, uint32_t number_threads )
{
    SyncList sync_list(number_threads);
    Filler filler(sync_list);
    std::list<Worker> worker_list;
    for( size_t i = 0; i < number_threads; ++i )
    {
        worker_list.emplace_back( Worker( i, sync_list ) );
    }

    // Start threads.
    auto filler_thread = std::async( std::launch::async, [&]{ filler.work(matrix, current_rank); } );
    std::vector< std::future<void> > worker_threads;
    for( auto& it: worker_list )
    {
        worker_threads.emplace_back( std::async( std::launch::async, [&]{ it.work(matrix); } ) );
    }

    // Wait to finish.
    for( auto& it: worker_threads )
    {
        it.get();
    }

    filler_thread.get();

    return current_rank;
}

