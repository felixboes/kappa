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
#include "parallelization.hpp"
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
    typedef typename MatrixType::MatrixEntryType MatrixEntryType;
    typedef typename MatrixType::DiagonalType DiagonalType;

    class Worker;

    /**
     *  Constructor.
    **/
    DiagonalizerField() : transp(false), def(0), rnk(0), num_working_threads(2), num_remaining_threads(0), current_rank(0) {}

    /**
     *  Diagonalizes a given matrix and gives access to the progress by writing the current rank to current_rank.
     *  If the number of working threads is greater than 1 and the number of remaining threads is greater than zero,
     *  we use the multithreaded version.
     *  @warning The list of rows we want to ommit has to be sorted.
    **/
    void operator() ( MatrixType& matrix );

    /**  @return defect of the matrix */
    uint32_t dfct() const;
    /**  @return defect of the matrix */
    HomologyField::KernT kern();
    /**  @return rank of the matrix */
    uint32_t rank() const;
    /**  @return rank of the matrix */
    HomologyField::TorsT tors();

//private:
    /**
     *  @return rank of matrix
     *  The matrix is diagonalized via Gauss to compute the number of linearly independant columns or rows.
     */
    uint32_t diag_field( MatrixType& matrix );

    /**
     *  @return rank of matrix and gives access to the progress by writing the current rank to current_rank.
     *  The matrix is diagonalized via Gauss to compute the number of linearly independant columns or rows.
     *  This version is parallelized and uses num_working_threads + num_remaining_threads many threads.
     *  For an explanation of the parallelization, see \c Worker.
     */
    uint32_t diag_field_parallelized( MatrixType& matrix );

    bool transp;    ///< True iff the transposed matrices are stored.
    uint32_t def;   ///< The defect of the matrix.
    uint32_t rnk;   ///< The rank of the matrix.
    uint32_t num_working_threads;       ///< The number of threads used in the process (diagoanlizing and collecting work).
    uint32_t num_remaining_threads;     ///< The number of threads collecting tasks.
    atomic_uint current_rank;           ///< The current state of the rank in the computation.
    
    std::list< size_t > ommit_rows;     ///< Rows we do not have to consider anymore.
    
    // Classes for paralellization:
    // Todo: give a detailed explanation.

    /**
     *  A JobQueue keeps track of the work that has to be done by all Workers.
     *  This work consists of:
     *  - A column and a row s.t. the corresponding matrix entry is non-zero.
     *    With this row, we want to erase all other entries in this column
     *    that do not belong to rows that already contribute to the rank.
     *  - A previously computed vector of all rows on which such
     *    row operations have to be performed, i.e. of the non-zero rows
     *    that do not yet contribute to the rank.
     *  - A vector of the zero rows of the matrix which do not yet
     *    contribute to the rank.
     *  When a new column is considered, the JobQueue can update this information.
     *  This assumes that the workers finish their work before.
     */
    class JobQueue
    {
    public:
        //! Collect initial work.
        JobQueue(
            MatrixType&                         matrix_init,
            typename MatrixType::DiagonalType&  diag,
            const std::list< size_t >&          ommit_these_rows,
            const uint32_t                      number_of_working_threads,
            const uint32_t                      number_of_remaining_threads,
            atomic_uint&                        current_rank );

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
        void get_work_chunk(uint32_t const thread_id, size_t & begin, size_t & end) const;
        void get_remaining_chunk(uint32_t const thread_id, size_t & begin, size_t & end) const;

        /*! Recomputes the chunk size for the working threads based on
        * num_working_threads and the current number of rows_to_work_at,
        * and the chunk size for the remaining threads based on num_remaining_threads
        * and the current number of remaining_rows.
        */
        void recompute_chunk_sizes();

        /*!
         * \brief Matrix to diagonalize.
         */
        MatrixType &              matrix;

        /*!
         * \brief Vector of rows r s.t. the entry (r, jobs.col) of the matrix is non-zero.
         * \note Rows that already contribute to the rank are not considered anymore.
         */
        std::vector<size_t>       rows_to_work_at;

        /*!
         * \brief Vector of rows r s.t. the entry (r, jobs.col) of the matrix is zero.
         * \note Rows that already contribute to the rank are not considered anymore.
         */
        std::vector<size_t>       remaining_rows;

        /*!
         * Number of rows in the vector rows_to_work_at
         * a single thread has to work at.
         */
        size_t                    work_chunk_size;
        size_t                    remaining_chunk_size;

        /*!
         * \brief Number of working_threads used.
         */
        size_t                    num_working_threads;

        /*!
         * \brief Number of remaining_threads used.
         */
        size_t                    num_remaining_threads;
        /*!
         * Row which is currently added to other rows for row operations.
         */
        size_t                    row_1;
        /*!
         * Current column.
         */
        size_t                    col;
    private:
        DiagonalType& diagonal;
    };

    /**
     * A Worker is supposed to perform row operations for a given column and to determine
     * the row operations for the next column.
     * More precisely:
     * At the beginning of the diagonalization of a matrix, we introduce
     * Workers and Threads of two different kinds explained later.
     * All Workers maintain a JobQueue of their current shared work.
     * At initialization, the JobQueue consists of
     *  - the 0th column of the matrix,
     *  - the first row row_1 for which the entry (col, row_1) is non-zero,
     *  - the vector of non-zero rows rows_to_work_at of the 0th column (except for row_1),
     *  - the vector of zero rows remaining_rows of the 0th column (except for row_1).
     * Now, we use all the Workers to either
     *  - perform the function work and thus perform row operations on a
     *    part of the rows_to_work_at, and thereby sort these rows into two new vectors,
     *    depending on whether their entry in the next column is zero or not,
     *  - sort the remaining_rows into two new vectors,
     *    depending on whether their entry in the next column is zero or not.
     * Afterwards, the function jobs.update_rank_and_work collects the rows
     * computed by the Workers, and they can continue their work until there is
     * no work left.
     * \note We assume that num_working_threads > 0 and num_remaining_threads > 0.
     */
    class Worker
    {
    public:
        Worker( const uint32_t identification, JobQueue & sl );

#ifdef BROKEN_VECTOR_IMPLEMENTATION
        Worker( Worker const & other);
        void operator=(Worker const & other);
#endif
        /*!
         * \brief Returns the id of this Worker, see \c id.
         */
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

        /*!
         * \brief Returns the new_rows_to_work_at.
         */
        std::vector<size_t> & get_new_rows_to_work_at();

        /*!
         * \brief Returns the new_remaining_rows.
         */
        std::vector<size_t> & get_new_remaining_rows();

    private:
        /*!
         * \brief id of this Worker
         * \note There are two kinds of ids: One for the working workers, one for the remaining workers.
         */
        uint32_t            id;

        /*!
         * \brief Contains all current job information of the workers.
         */
        JobQueue &          jobs;

        /*!
         * \brief Vector of rows r s.t. the entry (r, jobs.col + 1) of the matrix is non-zero.
         * This vector is built up by the worker.
         * \note Rows that already contribute to the rank are not considered anymore.
         */
        std::vector<size_t> new_rows_to_work_at;

        /*!
         * \brief Vector of rows r s.t. the entry (r, jobs.col + 1) of the matrix is zero.
         * This vector is built up by the worker.
         * \note Rows that already contribute to the rank are not considered anymore.
         */
        std::vector<size_t> new_remaining_rows;
    };
};

#include "diagonalizer_field.ipp"

#endif // DIAGONALIZER_FIELD_HPP
