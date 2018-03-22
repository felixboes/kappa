#include "permutation_manager.hpp"

PermutationManager::PermutationManager()
= default;

bool PermutationManager::is_trivial(const Transposition &t)
{
    return (t.first==t.second);
}

bool PermutationManager::have_disjoint_support(const Transposition &t1, const Transposition &t2)
{
    if(is_trivial(t1)  or  is_trivial(t2))
    {
        return true;
    }
    return !(t1.first==t2.first  or  t1.first==t2.second
             or t1.second==t2.first  or  t1.second==t2.second);
}

bool PermutationManager::are_symbolwise_disjoint(const Transposition &t1, const Transposition &t2)
{
    return !(t1.first==t2.first  or  t1.first==t2.second
             or t1.second==t2.first  or  t1.second==t2.second);
}

uint8_t PermutationManager::value_at(const Transposition &t, uint8_t n)
{
    uint8_t value_at_n = n;
    if(t.first==n)
    {
        value_at_n = t.second;
    }
    else if(t.second == n)
    {
        value_at_n = t.first;
    }
    return value_at_n;
}

uint8_t PermutationManager::value_at(const Norm2Permutation &t, uint8_t n)
{
    return value_at(t.first, value_at(t.second,n));
}

uint8_t PermutationManager::num_cycles(const Permutation& t)
{
    uint8_t p = (uint8_t) t.size() -1;
    uint8_t num_cycles = 0;
    std::vector<bool> visited( p+1, false );
    for( uint8_t i = 0; i <= p; )   // We iterate through all cycles and mark the visited symbols.
    {
        num_cycles += 1;
        visited.at(i) = true;
        uint8_t j = t.at(i);

        while( j != i and j <= p) // mark all symbols in this cycle
        {
            visited.at(j) = true;
            j = t.at(j);
        }
        if( j>p )
        {
            std::cerr << "Error in 'uint8_t PermutationManager::num_cycles(const Permutation& t) -> "
                      << "The Permutation is no permutation as some value is too large." << std::endl;
        }

        // find the next unvisited cycle
        for( ++i; i <= p and visited[i]; ++i )
        {
        }
    }
    return num_cycles;
}

std::vector < std::vector<uint8_t> >
PermutationManager::descending_cycle_decomposition_of_product(const Norm2Permutation &t1, const Norm2Permutation &t2)
{
    std::vector<uint8_t> occuring_symbols {t1.first.first, t1.first.second, t1.second.first, t1.second.second,
                                            t2.first.first, t2.first.second, t2.second.first, t2.second.second};
    std::sort(occuring_symbols.begin(), occuring_symbols.end());

    std::vector< std::vector<uint8_t> > cycle_decomp;
    std::map< uint8_t, bool > visited_symbols;

    for( int start_idx = occuring_symbols.size() -1; start_idx >= 0; ) // We iterate through all cycles
    {
        uint8_t start_symbol = occuring_symbols[start_idx];

        // determine the cycle of the start_symbol
        std::vector<uint8_t> cycle;
        uint8_t cur_symbol = start_symbol;

        do // find and store all symbols in this cycle, mark them as visited
        {
            cycle.push_back(cur_symbol);
            visited_symbols[cur_symbol] = true;
            cur_symbol = PermutationManager::value_at(t2, cur_symbol);
            cur_symbol = PermutationManager::value_at(t1, cur_symbol);
        }while(cur_symbol != start_symbol);

        // note that since the for-loop runs descendingly and occuring_symbols is sorted, the largest element of the
        // cycle is start_symbol, which is also cycle[0]
        if(cycle.size()>1)
        {
            cycle_decomp.push_back(cycle);
        }
        // find the next unvisited cycle
        for( --start_idx; start_idx >= 0; --start_idx )
        {
            if(visited_symbols.count(occuring_symbols[start_idx])==0) break;
        }
    }
    return cycle_decomp;
}

Permutation PermutationManager::inverse(const Permutation &t)
{
    Permutation t_inverse(t.size());
    for( uint8_t i = 0; i < t.size(); ++i )
    {
        t_inverse[ t.at(i) ] = i;
    }
    return t_inverse;
}

void PermutationManager::multiply(Permutation &sigma, const Transposition &t)
{
    std::swap( sigma[t.first], sigma[t.second] );
}

void PermutationManager::multiply(Permutation &sigma, const Norm2Permutation &t)
{
    multiply(sigma, t.first);
    multiply(sigma, t.second);
}

bool PermutationManager::has_norm_two(const Norm2Permutation &t)
{
    if (is_trivial(t.first) or is_trivial(t.second))
    {
        return false;
    }
    return !((t.first == t.second) or (t.first.first == t.second.second && t.first.second == t.second.first));
}

bool PermutationManager::write_in_norm2notation(Norm2Permutation &t)
{
    if( !has_norm_two(t) )
    {
        return false;
    }

    if(are_symbolwise_disjoint(t.first, t.second) )
    {
        if(t.first.first < t.first.second)
        {
            std::swap(t.first.first, t.first.second);
        }
        if(t.second.first < t.second.second)
        {
            std::swap(t.second.first, t.second.second);
        }
        if(t.first.first > t.second.first)
        {
            std::swap(t.first, t.second);
        }
        return true;
    }

    //remaining case: t is a threecycle; write t = (height(t), k, l)
    uint8_t height= std::max({t.first.first, t.first.second, t.second.first, t.second.second});
    uint8_t k = value_at(t, height);
    uint8_t l = value_at(t, k);
    t.first.first = l;
    t.first.second = k;
    t.second.first = height;
    t.second.second = l;
    if(k > l)
    {
        std::swap(t.first.first, t.first.second);
    }
    return true;
}

bool PermutationManager::is_in_norm2notation(const Norm2Permutation &t)
{
    return ( (t.first.first > t.first.second) and (t.second.first > t.second.second)
             and (t.second.first > t.first.first) );
}

void PermutationManager::substitute_index(Transposition &t, uint8_t idx, uint8_t replacement_idx)
{
    if(t.first == idx)
    {
        t.first = replacement_idx;
    }
    if(t.second == idx)
    {
        t.second = replacement_idx;
    }
}

void
PermutationManager::substitute_index(Norm2Permutation &t, uint8_t idx, uint8_t replacement_idx)
{
    substitute_index(t.first, idx, replacement_idx);
    substitute_index(t.second, idx, replacement_idx);
}

void PermutationManager::drop_fixed_point(Transposition &t, uint8_t fixed_pt)
{
    if(t.first > fixed_pt)
    {
        t.first = t.first -1;
    }
    if(t.second > fixed_pt)
    {
        t.second = t.second -1;
    }
}

void PermutationManager::drop_fixed_point(Norm2Permutation &t, uint8_t fixed_pt)
{
    drop_fixed_point(t.first, fixed_pt);
    drop_fixed_point(t.second, fixed_pt);
}

std::map<uint8_t, int8_t> PermutationManager::orientation_sign_for_ehrenfried(const Permutation &sigma)
{
    std::map< uint8_t, Permutation > cycle_decomp = sigma.cycle_decomposition();
    std::map< uint8_t, int8_t > sign;

    uint8_t i = 1; // counter of cycles
    for ( auto curr_cycle_it = cycle_decomp.begin(); curr_cycle_it != cycle_decomp.end(); ++curr_cycle_it )
    {
        Permutation cycle = curr_cycle_it->second;
        uint8_t min_symbol_of_cycle = curr_cycle_it->first;
        // if the cycle is a fixed point (a), we set sign(a) = 0 for the sake of completeness.
        if (cycle.is_fix_point(min_symbol_of_cycle))
        {
            sign[min_symbol_of_cycle] = 0;
            ++i;
            continue;
        }
        uint8_t second_min_symbol = 0;
        // determine the second min symbol of the cycle
        for (size_t m = 1; m < cycle.size(); ++m)
        {
            if (cycle.is_contained(m) && m != min_symbol_of_cycle)
            {
                second_min_symbol = m;
                break;
            }
        }

        // Let a_{j,1} be the first element in the j-th cycle in the decomposition; that is its minimum.
        // We have a_{1,1} < a_{2,1} < ... and a_{i,1} < second_min_symbol.
        // Find k such that a_{k,1} < second_min_symbol < a_{k+1,1}.
        uint8_t k = i;
        // Note that the iterator next_cycle_it always corresponds to the (k+1)-th cycle.
        // We found our position k iff the first time second_min_symbol < a_{k+1,1}.
        // In the other case we have again a_{k+1,1} < second_min_symbol. Induction.
        for ( auto next_cycle_it = (std::next(curr_cycle_it)); next_cycle_it != cycle_decomp.end(); ++next_cycle_it )
        {
            uint8_t next_min_symbol = next_cycle_it->first;
            if ( second_min_symbol < next_min_symbol )
            {
                if ( ((k - i) % 2) == 0 )
                {
                    sign[min_symbol_of_cycle] = 1;
                }
                else
                {
                    sign[min_symbol_of_cycle] = -1;
                }
                break;
            }
            ++k;
        }
        // If second_min_symbol > a_{m, 1}, we still need to determine the sign of min_symbol_of_cycle.
        if (k == cycle_decomp.size())
        {
            if ( ((k - i) % 2) == 0 )
            {
                sign[min_symbol_of_cycle] = 1;
            }
            else
            {
                sign[min_symbol_of_cycle] = -1;
            }
        }

        // for all other symbols of the cycle, we set the sign to 1
        for (size_t c = 0; c < cycle.size(); ++c)
        {
            if (cycle.is_contained(c) and c != min_symbol_of_cycle)
            {
                sign[c] = 1;
            }
        }
        ++i;
    }
    return sign;
}

