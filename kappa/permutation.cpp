#include "permutation.hpp"

Permutation::Permutation()
    : data()
{
    // intentionally do nothing
}

Permutation::Permutation(const uint8_t size)
    : data(size, size)
{
    // intentionally do nothing
}

Permutation::Permutation(const Permutation & other)
    :
      data(other())
{
    // intentionally do nothing
}

std::vector<uint8_t> Permutation::operator()() const
{
    return data;
}

uint8_t & Permutation::operator[](const uint8_t i)
{
    return data[i];
}

uint8_t const & Permutation::at(const uint8_t i) const
{
    return data.at(i);
}

uint8_t Permutation::size() const
{
    return data.size();
}

bool Permutation::operator==(const Permutation &t) const {
    return (data==t());
}

std::map< uint8_t, Permutation > Permutation::cycle_decomposition () const
{
    const uint8_t p = this->size() - 1;
    std::map< uint8_t, Permutation > cycle_decomp;
    std::vector< bool > visited( p+1, false );
    for( uint8_t i = 0; i <= p; ) // We iterate through all cycles and mark the used symbols.
    {
        // determine the cycle of i.
        Permutation cycle(p+1);

        uint8_t prev;     // previous symbol
        uint8_t cur = i; // current symbol

        do // mark all symbols in this cycle
        {
            visited[cur] = true;
            prev = cur;
            cur = this->at(prev);
            cycle[prev] = cur;
        }while( cur != i );
        // note that since the for-loop runs ascendingly, the smallest element of the cycle
        // is i
        cycle_decomp[i] = cycle;
        // find the next unvisited cycle
        for( ++i; i <= p && visited[i]; ++i )
        {
        }
    }
    return cycle_decomp;
}

bool Permutation::is_contained(const uint8_t i) const
{
    return (data[i] != data.size());
}

bool Permutation::is_fix_point(const uint8_t i) const
{
    return (data[i] == i);
}

std::ostream& operator<< (std::ostream& stream, const Permutation& permutation)
{
    stream << "permutation: " << std::endl;
    for (uint8_t i = 0; i < permutation.size(); ++i)
    {
        if (permutation.is_contained(i))
        {
            stream << (size_t) i << " maps to " << (size_t) permutation.at(i) << std::endl;
        }
    }
    return stream;
}

Permutation Permutation::long_cycle(uint8_t p)
{
    Permutation sigma(p+1);
    for(uint8_t k = 0; k < p; ++k)
    {
        sigma[k] = k+1;
    }
    sigma[p] = 0;
    return sigma;
}

Permutation Permutation::long_cycle_inv(uint8_t p)
{
    Permutation sigma(p+1);
    for(uint8_t k = 1; k <= p; ++k)
    {
        sigma[k] = k-1;
    }
    sigma[0] = p;
    return sigma;
}
