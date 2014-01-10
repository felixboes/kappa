#include <chrono>
#include <future>
#include <iostream>
#include <list>
#include <mutex>

struct RowOpParam
{
    RowOpParam(size_t r1 = 0, size_t r2 = 0, size_t c = 0) : row1(r1), row2(r2), col(c) {}
    size_t row1;
    size_t row2;
    size_t col;
};

class SyncList
{
public:
    SyncList() : genug_fuer_heute(false) {}
    void refill(size_t col);
    bool get(RowOpParam &);

    void ja_genug_fuer_heute() { std::lock_guard<std::mutex> lock(mtx); genug_fuer_heute = true; cond.notify_all(); }
    bool keine_arbeit_mehr() { std::lock_guard<std::mutex> lock(mtx); return genug_fuer_heute; }
private:
    std::list<RowOpParam> row_op_list;
    std::mutex mtx;
    std::condition_variable cond;
    bool genug_fuer_heute;
};

void SyncList::refill(size_t col)
{
    std::lock_guard<std::mutex> lock(mtx);
    for( size_t i = 1; i < 10; ++i )
    {
        row_op_list.push_back( RowOpParam(i, 2*i, col) );
    }
    cond.notify_all();
}

bool SyncList::get( RowOpParam & rop )
{
    std::unique_lock<std::mutex> lock(mtx);
    cond.wait( lock, [this]{ return ! ( row_op_list.empty() && !genug_fuer_heute ) ; } );
    if( !row_op_list.empty() )
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




class Worker
{
public:
    Worker( uint32_t nummer, SyncList & wl ) : id(nummer), work_list(wl) {}
    void work();
private:
    uint32_t id;
    SyncList& work_list;
};

void Worker::work()
{
    RowOpParam rop;
    
    while( work_list.get(rop) == true )
    {
        std::cout << id << ": " << rop.row1 << " " << rop.row2 << " " << rop.col << std::endl;
    }
}





int main()
{
    SyncList arbeitsliste;

    auto arbeitgeber = std::async( std::launch::async, [&]()
    {
        for( int i = 0; i < 2; ++i)
        {
            std::cout << "Fill" << std::endl;
            arbeitsliste.refill(i);
            std::this_thread::sleep_for( std::chrono::seconds(1) );
        }
        std::cout << "Filled." << std::endl;
        arbeitsliste.ja_genug_fuer_heute();
    } );
    
    Worker arbeiter1 (1, arbeitsliste);
    Worker arbeiter2 (2, arbeitsliste);
    
    auto arbeit1 = std::async( std::launch::async, [&]{ arbeiter1.work(); } );
    auto arbeit2 = std::async( std::launch::async, [&]{ arbeiter2.work(); } );
    
    arbeit1.get();
    arbeit2.get();
    arbeitgeber.get();
    
    return 0;
}
