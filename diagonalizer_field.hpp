#ifndef DIAGONALIZER_FIELD_HPP
#define DIAGONALIZER_FIELD_HPP

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
 *  coefficients in \f$ \mathbb{F}_p \f$).
 */

template < class CoefficientT >
class DiagonalizerField
{
public:
    typedef CoefficientT CoefficientType;
    typedef MatrixField<CoefficientType> MatrixType;
    
    /**  constructor */
    DiagonalizerField() {}
    
    /**
    *   Diagonalize a given matrix.
    **/
    void operator() ( MatrixType &matrix, uint32_t number_threads=0 );
    void operator() ( MatrixType &matrix, atomic_uint & current_rank, uint32_t number_threads=0 );
    
    /**
     *  Diagonalize matrix and apply the base change to post_matrix.
     *  One should have the following picture in mind.
     *  \f[
     *      C_{n-2} \xleftarrow{\partial_{postmatrix}} C_{n-1} \xleftarrow{\partial_{matrix}} C_n
     *  \f]
     */
    void operator() ( MatrixType &post_matrix, MatrixType &matrix );
    
    /**
     *  Diagonalize matrix and apply the base change to pre_matrix and post_matrix.
     *  One should have the following picture in mind.
     *  \f[
     *      C_{n-2} \xleftarrow{\partial_{postmatrix}} C_{n-1} \xleftarrow{\partial_{matrix}} C_n \xleftarrow{\partial{prematrix}} C_{n+1}
     *  \f]
     */
    void operator() ( MatrixType &post_matrix, MatrixType &matrix, MatrixType &pre_matrix );
    
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
    uint32_t diag_field(MatrixType& matrix, atomic_uint & current_rank);
    
    uint32_t diag_field_parallelized(MatrixType& matrix, atomic_uint & current_rank, uint32_t number_threads);
    
    MatrixType in;
    MatrixType out;
    uint32_t def;
    uint32_t rnk;
     
    // Classes for paralellization:
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
