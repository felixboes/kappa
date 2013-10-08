#ifndef HOMOLOGY_FIELD_HPP
#define HOMOLOGY_FIELD_HPP

#include <cinttypes>
#include <iomanip>
#include <iostream>
#include <map>

class HomologyField
{
public:
    typedef int64_t KernT;
    typedef int64_t TorsT;
    
    HomologyField();
    HomologyField( int32_t, KernT, TorsT );
    
    friend std::ostream& operator<< (std::ostream& stream, const HomologyField& homol);
private:
    std::map< int32_t, int64_t > dimension;
};

std::ostream& operator<< (std::ostream& stream, const HomologyField& homol);

#endif // HOMOLOGY_FIELD_HPP
