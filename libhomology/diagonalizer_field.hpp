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
     *   A RowOpParam can be seen as a job that has to be processed by one of the worker threads.
     */
    struct RowOpParam
    {
        RowOpParam(size_t r1 = 0, size_t r2 = 0, size_t c = 0) : row1(r1), row2(r2), col(c) {}
        size_t row1;
        size_t row2;
        size_t col;
    };

    /**
     *  A JobQueue keeps track of the work that has to be done.
     *  It can fill itself with jobs i.e. RowOpParam. Furthermore it can tell a
     *  thread which chunk of jobs are to be performed by this thread.
     */
    class JobQueue
    {
    public:
        //! Basic constructors
        JobQueue(MatrixType & matrix_init, uint32_t number_of_working_threads );


#ifdef BROKEN_VECTOR_IMPLEMENTATION
        JobQueue(JobQueue const & other);
        JobQueue & operator=(JobQueue const & other);
#endif

        //! Access operation with a given id
        RowOpParam const & get_operation(size_t op_id) const;

        //! Get the chunk [begin, end) of operations dedicated to the given thread
        void get_chunk(uint32_t const thread_id, size_t & begin, size_t & end) const;

        /*! \brief Computes new operations and returns whether there is new work
         *
         *  The former operations will be purged.
        **/
        bool compute_operations(atomic_uint & cur_rank);

    private:
        //! Recomputes the chunk size based on the current number of operations
        void recompute_chunk_size();

        MatrixType &              matrix;
        std::vector<RowOpParam>   operations;

        size_t                    cur_chunk_size;
        size_t                    number_of_threads;

        std::vector<size_t>       rows_to_check; //!< Auxiliary for recompute_chunk_size
        size_t                    col;
    };

    /**
     *  A Worker waits for row operation jobs and performs those.
     */
    class Worker
    {
    public:
        Worker( uint32_t identification, JobQueue & sl );

#ifdef BROKEN_VECTOR_IMPLEMENTATION
        Worker( Worker const & other);
        Worker & operator=(Worker const & other);
#endif


        //! Perform the thread dedicated work. (As determined by the JobQueue.)
        void work(MatrixType & matrix);

    private:
        uint32_t   id;
        JobQueue & jobs;
    };
};

#include "diagonalizer_field.ipp"

#endif // DIAGONALIZER_FIELD_HPP
