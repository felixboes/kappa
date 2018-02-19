#include "homology.hpp"


int main()
{
    typedef DiagonalizerField<Zm> DFZm;
    DFZm::SyncList arbeitsliste;

    auto arbeitgeber = std::async( std::launch::async, [&]()
    {
        for( int i = 0; i < 2; ++i)
        {
            std::cout << "Fill" << std::endl;
            arbeitsliste.refill(i);
            std::this_thread::sleep_for( std::chrono::seconds(1) );
        }
        std::cout << "Filled." << std::endl;
        arbeitsliste.all_work_done();
    } );
    
    DFZm::Worker arbeiter1 (1, arbeitsliste);
    DFZm::Worker arbeiter2 (2, arbeitsliste);
    
    auto arbeit1 = std::async( std::launch::async, [&]{ arbeiter1.work(); } );
    auto arbeit2 = std::async( std::launch::async, [&]{ arbeiter2.work(); } );
    
    arbeit1.get();
    arbeit2.get();
    arbeitgeber.get();
    
    return 0;
}
