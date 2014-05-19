#ifndef DIAGONALIZER_FIELD_HPP
#define DIAGONALIZER_FIELD_HPP

// Description:
//
// This header defines a function object that diagonalizes matrices with field coefficients 'CoefficientT'.
// Here we offer single- and multithreaded versions.

#include <chrono>
#include <cinttypes>
#include <functional>
#include <future>
#include <list>
#include <mutex>
#include <thread>
#include <vector>

#include "homology_field.hpp"
#include "matrix_field.hpp"
#include "thread.hpp"

/**
 *  In order to compute the homology of the chain complex
 *  \f[
 *          C_{n-1} \xleftarrow{\partial_{out}} C_n \xleftarrow{\partial_{in}} C_{n+1}
 *  \f]
 *  at position \f$n\f$, it suffices to compute the dimension of the kernel and the image of the respective matrices.
 *  Therefore we need to diagonalize the matrices representing the two differentials (here only for matrices with
 *  coefficients in \f$\mathbb Q\f$ and \f$ \mathbb{F}_p \f$).
 */
template < class MatrixT >
class DiagonalizerField
{
public:
    typedef MatrixT MatrixType;    ///< We use this typedef to grant access the matrix type from other classes.
    class Worker;

    /**
     *  Constructor.
     */
    DiagonalizerField() {}

    /**
    *   Diagonalizes a given matrix.
    *   If the number of threads is given and greater then 1, we use the multithreaded version.
    **/
    void operator() ( MatrixType &matrix, uint32_t number_threads=0 );

    /**
    *   Diagonalizes a given matrix and gives access to the progress by writing the current rank to current_rank.
    *   If the number of threads is given and greater then 1, we use the multithreaded version.
    **/
    void operator() ( MatrixType &matrix, atomic_uint & current_rank, uint32_t number_threads=0 );

    /**  @return defect of the matrix */
    uint32_t dfct();
    /**  @return defect of the matrix */
    HomologyField::KernT kern();
    /**  @return rank of the matrix */
    uint32_t rank();
    /**  @return rank of the matrix */
    HomologyField::TorsT tors();

//private:
    /**
     *  @return rank of matrix
     *  The matrix is diagonalized via Gauss to compute the number of linearly independant columns or rows.
     */
    uint32_t diag_field(MatrixType& matrix);

    /**
     *  @return rank of matrix and gives access to the progress by writing the current rank to current_rank.
     *  The matrix is diagonalized via Gauss to compute the number of linearly independant columns or rows.
     */
    uint32_t diag_field(MatrixType& matrix, atomic_uint & current_rank);

    /**
     *  @return rank of matrix and gives access to the progress by writing the current rank to current_rank.
     *  The matrix is diagonalized via Gauss to compute the number of linearly independant columns or rows.
     *  This version is parallelized and uses number_threads many threads.
     */
    uint32_t diag_field_parallelized(MatrixType& matrix, atomic_uint & current_rank, uint32_t number_threads);

    uint32_t def;   ///< The defect of the matrix.
    uint32_t rnk;   ///< The rank of the matrix.

    // Classes for paralellization:
    // Todo: give a detailed explanation.

    /**
     *  A JobQueue keeps track of the work that has to be done.
     *  This work consists of:
     *  - A column and a row s.t. the corresponding matrix entry is non-zero.
     *    With this row, we want to erase all other entries in this column
     *    that do not belong to rows that already contribute to the rank.
     *  - A previously computed vector of all rows on which such
     *    row operations have to be performed.
     *  - A vector of the remaining rows of the matrix which do not yet
     *    contribute to the rank.
     *  When a new column is considered, the JobQueue can update this information.
     *  This assumes that the workers finish their work before.
     */
    class JobQueue
    {
    public:
        //! Collect initial work.
        JobQueue(MatrixType & matrix_init, uint32_t number_of_working_threads, atomic_uint & current_rank );

#ifdef BROKEN_VECTOR_IMPLEMENTATION
        JobQueue(JobQueue const & other);
        void operator=(JobQueue const & other);
#endif

        /**
         * This function is supposed to be called by update_rank_and_work.
         * It fetches the all work for the next column computed by the threads
         * and unites it to become the current work.
         * ToDo move into c file
         */
        void add_up_rows
            (std::vector<Worker> & workers,
             std::vector<size_t> & former_rows_to_work_at,
             std::vector<size_t> & former_remaining_rows);

        /**
         * This function is supposed to be called after all threads finish their work.
         * It fetches the all work for the next column computed by the threads
         * and unites it to become the current work.
         * It also sets the new row_1, rank and chunk size accordingly.
         * @return True iff there is work left.
         */
        bool update_rank_and_work(atomic_uint & current_rank, std::vector<Worker> & workers);

        //! Get the chunk [begin, end) of rows_to_work_at dedicated to the given thread.
        void get_chunk(uint32_t const thread_id, size_t & begin, size_t & end) const;

        //! Recomputes the chunk size based on the current number of rows_to_work_at.
        void recompute_chunk_size();

        MatrixType &              matrix;
        /*!
         * Rows of the matrix on which row operations w.r.t. col and row_1
         * have to be performed.
         */
        std::vector<size_t>       rows_to_work_at;
        /*!
         * All rows of the matrix that do not yet contribute to the rank
         * and that are not within the vector rows_to_work_at.
         */

        std::vector<size_t>       remaining_rows;

        /*!
         * Number of rows in the vector rows_to_work_at
         * a single thread has to work at.
         */
        size_t                    cur_chunk_size;
        /*!
         * Total number of threads used for parallelization.
         * One thread takes care of the remaining_rows while all the others
         * perform row operations on the rows_to_work_at.
         * TODO Check whether this can be improved!
         */
        size_t                    number_of_threads;
        /*!
         * Row which is currently added to other rows for row operations.
         */
        size_t                    row_1;
        /*!
         * Current column.
         */
        size_t                    col;
    };

    /**
     * A Worker is supposed to perform row operations for a given column and to determine
     * the row operations for the next column.
     * More precisely:
     * All Workers maintain a JobQueue of their current shared work.
     * At initialization, the JobQueue consists of
     *  - the 0th column of the matrix,
     *  - the first row for which the entry (col, row_1) is non-zero,
     *  - the set of non-zero rows rows_to_work_at of the 0th column (except for row_1),
     *  - the set of remaining rows remaining_rows of the 0th column (except for row_1).
     * Now, a Worker can either
     *  - perform the function work and thus performs row operations on a
     *    part of the rows_to_work_at, and thereby sort these rows into two new vectors,
     *    depending on whether their entry in the next column is zero or not,
     *  - sort the remaining_rows into two new vectors,
     *    depending on whether their entry in the next column is zero or not.
     * The Workers must be distributed s.t. all rows are considered.
     * Afterwards, the function jobs.update_rank_and_work collects the rows
     * computed by the Workers, and they can continue their work until there is
     * no work left.
     */
    class Worker
    {
    public:
        Worker( uint32_t identification, JobQueue & sl );

#ifdef BROKEN_VECTOR_IMPLEMENTATION
        Worker( Worker const & other);
        void operator=(Worker const & other);
#endif
        uint32_t get_id()
        {
            return id;
        }

        /*!
         * Perform row operations on the corresponding chunk of the
         * rows_to_work_at, and thereby sort these rows into two new vectors,
         * depending on whether their entry in the next column is zero or not.
         */
        void work(MatrixType & matrix);

        /*!
         * Sort the remaining_rows into two new vectors,
         * depending on whether their entry in the next column is zero or not.
         */
        void collect_remaining_work(MatrixType & matrix);

        std::vector<size_t> & get_new_rows_to_work_at();

        std::vector<size_t> & get_new_remaining_rows();

    private:
        uint32_t            id;
        JobQueue &          jobs;
        std::vector<size_t> new_rows_to_work_at;
        std::vector<size_t> new_remaining_rows;
    };
};

#include "diagonalizer_field.ipp"

#endif // DIAGONALIZER_FIELD_HPP
