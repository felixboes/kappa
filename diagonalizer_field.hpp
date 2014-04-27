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
     *  A SyncList keeps track of the work that has to be done.
     *  It can be filled with jobs i.e. RowOpParam.
     */
    class SyncList
    {
    public:
        SyncList(uint32_t number_of_working_threads ) : num_busy_workers(number_of_working_threads), work_done(false) {}
        
        /**
         *  blocks till queue is empty.
         */
        void wait_till_empty();
        
        /**
         *  puts a new job into the queue.
         */ 
        void put(RowOpParam new_job);
        
        /**
         *  get stores a new job in new_job or blocks if the queue is empty.
         *  get returns true iff there was a new job.
         */
        bool get(RowOpParam & new_job);
        
        /**
         *  Tell the waiting workers that there is new work to do, if the number of workers is n.
         */
        void wait_for_workers();
        
        /**
         *  Tell the SyncList that all work is done.
         */
        void all_work_done();
        
        /**
         * @returns true iff there is no more work to do.
         */
        bool no_work_left();
    private:
        std::list<RowOpParam> row_op_list;
        std::mutex mtx;
        std::condition_variable cond_wait_for_filler;
        std::condition_variable cond_wait_for_worker;
        uint32_t num_busy_workers;
        bool work_done;
    };
    
    /**
     *  A Filler refills the SyncList.
     */
    class Filler
    {
    public:
        Filler( SyncList & sl );
        
        /**
         *  Refill the SyncList with new work evey time a column is fully processed.
         */
        void work(MatrixType& matrix, atomic_uint & current_rank);
    private:
        SyncList & sync_list;
    };
    
    /**
     *  A Worker waits for row operation jobs and performs those.
     */
    class Worker
    {
    public:
        Worker( uint32_t identification, SyncList & sl );
        
        /**
         *  Wait and perform row operations till all work is done.
         */ 
        void work(MatrixType& matrix);
    private:
        uint32_t id;
        SyncList& sync_list;
    };
};

#include "diagonalizer_field.ipp"

#endif // DIAGONALIZER_FIELD_HPP
