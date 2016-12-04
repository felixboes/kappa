#include "tuple.hpp"

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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Tuple::radial = false;
uint32_t Tuple::min_symbol = 1;
uint32_t Tuple::min_boundary_offset = 1;
uint32_t Tuple::max_boundary_offset = 1;

Tuple :: Tuple() :
    p(0),
    rep()
{
}

Tuple :: Tuple(const size_t h) :
    p(0),
    rep( h, Transposition(0, 0) )
{
}

Tuple :: Tuple(const uint32_t symbols, const size_t h) :
    p(symbols),
    rep( h, Transposition(0, 0) )
{
}

void Tuple::parallel_case()
{
    radial = false;
    min_symbol = 1;
    min_boundary_offset = 1;
    max_boundary_offset = 1;
}

void Tuple::radial_case()
{
    radial = true;
    min_symbol = 0;
    min_boundary_offset = 0;
    max_boundary_offset = 0;
}

bool Tuple::get_radial()
{
    return radial;
}

uint32_t Tuple::get_min_symbol()
{
    return min_symbol;
}

uint32_t Tuple::get_min_boundary_offset()
{
    return min_boundary_offset;
}

uint32_t Tuple::get_max_boundary_offset()
{
    return max_boundary_offset;
}

Transposition & Tuple :: at(const size_t n)
{
    return rep[n-1];
}

const Transposition & Tuple :: at(const size_t n) const
{
    return rep[n-1];
}

Transposition& Tuple :: operator[](const size_t n)
{
    return at(n);
}

int32_t Tuple :: norm() const
{
    if( this->operator bool() == true )
    {
        return rep.size();
    }
    else // The Tuple is not valid.
    {
        return 0;
    }
}

bool Tuple :: operator==(const Tuple& t) const
{
    // Observe how two vectors are compared
    // Operations == and != are performed by first comparing sizes, and if they match, the elements are compared sequentially
    // using algorithm equal, which stops at the first mismatch.
    // Source: http://www.cplusplus.com/reference/vector/vector/operators/
    return (this->rep == t.rep);
}

bool Tuple :: operator!=(const Tuple& t) const
{
    // Observe how two vectors are compared
    // Operations == and != are performed by first comparing sizes, and if they match, the elements are compared sequentially
    // using algorithm equal, which stops at the first mismatch.
    // Source: http://www.cplusplus.com/reference/vector/vector/operators/
    return (this->rep != t.rep);
}

Tuple :: operator bool() const
{
    if (rep.size() == 0)
    {
        return false;
    }
    return true;
}

bool Tuple :: has_correct_num_cycles(const size_t m) const
{
    // min_symbol = 1 iff radial = false iff num_components = m + 'exactly one boundary component of the surface'
    return num_cycles() == m + min_symbol;
}

bool Tuple :: is_multiple_of_a() const
{
    if( this->operator bool() == true )
    {
        if( this->at( this->rep.size() ) == Transposition( this->p, this->p-1 ) )
        {
            const size_t rep_size_min_1 = rep.size() - 1;
            const size_t p_min_2 = this -> p - 2;
            for( size_t i = 0; i < rep_size_min_1; ++i )
            {
                if( rep.at(i).first > p_min_2 || rep.at(i).second > p_min_2 )
                {
                    return false;
                }
            }
            return true;
        }
    }

    return false;
}

uint32_t Tuple :: num_cycles() const
{
    // Since the t_i are transpositions, one can instead count the number of cycles of (p p-1 ... 1) t_1 ... t_h.
    uint32_t num_cycles = 0;
    Permutation sigma_inv = long_cycle_inv();

     // multiply with t_1, ..., t_h
    for( int32_t i = 1; i <= norm(); i++ )
    {
        uint8_t k = at(i).first;
        uint8_t l = at(i).second;
        std::swap( sigma_inv[ k ], sigma_inv[ l ] );
    }

    // count the cycles
    std::vector<bool> visited( p+1, false );    // visited[0] is not used
    for( uint8_t i = min_symbol; i <= p; ) // We iterate through all cycles and mark the used symbols.
    {
        // consider the next cycle
        num_cycles += 1;
        visited[i] = true;
        uint8_t j = sigma_inv[i];
        
        while( j != i ) // mark all symbols in this cycle
        {
            visited[j] = true;
            j = sigma_inv[j];
        }

        // find the next unvisited cycle
        for( ++i; i <= p && visited[i]; ++i )
        {
        }
    }
    return num_cycles;
}

Tuple::ConnectedComponents Tuple::connected_components() const
{
    // TODO Adapt for radial case
    // Compare http://www.boost.org/doc/libs/1_49_0/libs/graph/example/connected_components.cpp
    typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS> Graph;

    // Build graph.
    Graph G;
    for( auto edge : rep )
    {
        boost::add_edge(edge.first, edge.second, G);
    }

    // Let boost compute the connected compontents.
    ConnectedComponents components(p+1);
    boost::connected_components(G, &components[0]);

    return components;
}

/*
std::vector< size_t > Tuple::shuffle_positions() const
{
    std::vector< size_t > slits;
    const size_t h = norm();
    size_t k = 2*h;
    slits.push_back(k);
    
    while( k != 0 )
    {
        const size_t j = (k+1)/2;
        const size_t a = ( 2*j - k == 1 ? at(j).first : at(j).second );
        
        std::cout << j << " -> " << a << std::endl;
        
        // search for slits on the same hight but below a (i.e. longer slits).
        for( size_t l = j-1; l >= 1; --l )
        {
            if( a == at(l).first )
            {
                k = 2*l;
                // [n3290: 6.1/1]: [..] The scope of a label is the function in which it appears. [..]
                // Therefore we can use the name 'end_outer_loop'.
                goto end_outer_loop;
            }
            else if( a == at(l).second )
            {
                k = 2*l-1;
                // [n3290: 6.1/1]: [..] The scope of a label is the function in which it appears. [..]
                // Therefore we can use the name 'end_outer_loop'.
                goto end_outer_loop;
            }
        }
        
        // if a is the longest slit at the first line, we are done.
        if( a == 1 )
        {
            k = 0;
            // [n3290: 6.1/1]: [..] The scope of a label is the function in which it appears. [..]
            // Therefore we can use the name 'end_outer_loop'.
            goto end_outer_loop;
        }
        
        // find the first slit below a.
        for( size_t l = h; l >= 1; --l )
        {
            if( a-1 == at(l).first )
            {
                k = 2*l;
                // [n3290: 6.1/1]: [..] The scope of a label is the function in which it appears. [..]
                // Therefore we can use the name 'end_outer_loop'.
                goto end_outer_loop;
            }
            else if( a-1 == at(l).second )
            {
                k = 2*l-1;
                // [n3290: 6.1/1]: [..] The scope of a label is the function in which it appears. [..]
                // Therefore we can use the name 'end_outer_loop'.
                goto end_outer_loop;
            }
        }
        
        end_outer_loop:
        slits.push_back(k);
    }
    return slits;
}

std::vector< size_t > Tuple::slits() const
{
    std::vector< size_t > the_slits;
    size_t h = norm();
    
    for( size_t i = 1; i <= p; ++i )
    {
        for( size_t j = 1; j <= h; ++j )
        {
            if( at(j).first == i )
            {
                the_slits.push_back( 2*j );
            }
            else if ( at(j).second == i )
            {
                the_slits.push_back( 2*j - 1 );
            }
        }
    }
    
    return the_slits;
}

Tuple Tuple::Q_term( const std::vector< size_t >& shuffle_slit_conf, const size_t num_shuffle_pos, const size_t fuse_pos ) const
{
    const size_t h = norm();
    Tuple res(2*p, 2*norm()); // We reduce p later.
    const auto the_slits = slits();
    const auto the_shuffles = shuffle_positions();
    
    for( size_t j = 0; j < 2*h; ++j )
    {
                                                                         
    }
    
    return res;
}
*/

int32_t Tuple::num_clusters() const
{
    typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS> Graph;

    // Build graph.
    Graph G;
    for( auto edge : rep )
    {
        boost::add_edge(edge.first, edge.second, G);
    }

    // Let boost compute the connected compontents.
    ConnectedComponents components(p+1);
    // returns the number of components.
    // Here 0 is a valid vertex, so we must shrink the number of components if radial = false.
    return boost::connected_components(G, &components[0]) - min_symbol;
}

bool Tuple :: monotone() const
{
    // A tuple is monotone iff the sequence of all at(i) is monotone.
    for( int32_t i = 1; i <= norm() - 1; i++ )
    {
        if( at(i+1).first < at(i).first )
        {
            return false;
        }
    }
    return true;
}

bool Tuple :: f(const uint32_t i)
{
    // Denote g_{i+1} | g_i by (ab)(cd).
    uint8_t a = at(i+1).first;
    uint8_t b = at(i+1).second;
    uint8_t c = at(i).first;
    uint8_t d = at(i).second;

    if( a == c ) // (a )(a )
    {
        if( b == d ) // (ab)(ab) = id and (a*)(a*) = 0
        {
            rep = std::vector<Transposition>();
            return false;
        }
        else if( a == b ) // (a*)(ad) = (ad*) -> (d*)(ad)
        {
            at(i+1) = Transposition(d,d);
            return true;
        }
        else if( a == d ) // (ab)(a*) = (ab*) -> (b*)(ab)
        {
            at(i+1) = Transposition(b,b);
            at(i)   = Transposition(a,b);
            return true;
        }
        else // (ab)(ad) = (adb) = (bd)(ab) or (db)(ab)
        {
            at(i+1) = Transposition(std::max(b,d), std::min(b,d));
            at(i)   = Transposition(a,b);
            return true;
        }
    }
    else if( a == b ) // (a*)(  )
    {
        if( a == d ) // (a*)(ca) = (ca*) -> (a*)(ca)
        {
            return true;
        }
        else // (a*)(c*) or (a*)(cd) and the two transpositions commute
        {
            if( a < c )
            {
                return true;
            }
            else
            {
                std::swap( at(i+1), at(i));
                return true;
            }
        }
    }
    // Case: a > b
    else if( b == c ) // (ab)(b )
    {
        if( c == d ) // (ab)(b*) = (ab*) -> (b*)(ab)
        {
            std::swap( at(i+1), at(i) );
            return true;
        }
        else // (ab)(bd) = (abd) and (abd)' = (ad) and (abd)(ad) = (bd)
        {
            at(i+1) = Transposition(b,d);
            at(i)   = Transposition(a,d);
            return true;
        }
        // These are all cases since a > b >= d.
    }
    else if( a == d ) // (ab)(ca) = (acb) and c > a -> (ab)(ca)
    {
        return true;
    }
    else if( b == d ) // (ab)(cb) = (abc) -> a < c: (ab)(cb) oder a > c:(cb)(ac)
    {
        if( a < c )
        {
            return true;
        }
        else    // (abc)' = (ac) and (abc)(ac) = (cb)
        {
            at(i+1) = Transposition(c,b);
            at(i)   = Transposition(a,c);
            return true;
        }
    }
    else // (ab)(c ) with d not a and not b, hence disjoint transpositions
    {
        if( a < c )
        {
            return true;
        }
        else
        {
            std::swap( at(i+1), at(i) );
            return true;
        }
    }
    std::cerr << "Error in 'bool Tuple::f(uint32_t i)' -> Reached impossible case." << std::endl;
    return false;
}
                 
bool Tuple :: phi( const uint32_t q, const uint32_t i )
{
    if( i == 0 || i > q )
    {
        std::cerr << "Error in 'bool Tuple::phi( uint32_t q, uint32_t i)' -> q=" << q << " and i=" << i << std::endl;
        return false;
    }

    for( uint32_t j = q-1; j >= i; j-- ) // The loop terminates due to i > 0.
    {
        if( f(j) == false ) // The norm of the product falls.
        {
            return false;
        }
    }
    return true;
}

Tuple Tuple :: d_hor( const uint8_t k ) const
{
    Tuple boundary = *this;
    
    // start with sigma_0.
    Permutation sigma = long_cycle();
    Permutation sigma_inv = long_cycle_inv();
    
    for(uint8_t q = 1; q <= boundary.norm(); ++q)
    {
        // Write tau_q = (a,b)
        auto a = boundary[q].first;
        auto b = boundary[q].second;
        
        //Write (k, sigma_{q-1}(k)) = (k,l)
        auto l = sigma[k];
        
         // Compute tau':
        // Most of the time the transpositions are disjoint hence (a,b)(k,l) = (k,l)(a,b) and
        // the left transposition will be killed by D_k
        if( k != a && k != b && l != a && l != b )
        {
        }
        // The degenerate case:
        // k and l are both part of tau_q.
        else if( a == std::max(k, l) and ( b == std::min(k, l) or k == l))
        {
            return Tuple();
        }
        // The non degenerate case:
        else
        {   
            // Compute Z = (a,b)(k,l)
            // Z(k) = l iff l != a and l != b hence k == a or k == b
            // In this case: (Z(k),k)Z = (k,l)(a,b)(k,l) = (c,l) with c != k
            if( l != a && l != b)
            {
                if(a != k)
                {
                    boundary[q].first  = std::max(a,l);
                    boundary[q].second = std::min(a,l);
                }
                else
                {
                    boundary[q].first  = std::max(b,l);
                    boundary[q].second = std::min(b,l);
                }
            }
            // Z(k) != l iff l = a or l = b hence k != a and k != b
            // In this case: (Z(k),k)Z = (c,k)(a,b)(k,l) = (c,l) with c != l,
            // but we see that this is just (a,b):
            // if l != a, then k,l != a and we map
            // a to a to b != c since otherwise a would map to k
            // therefore a = c and l = b.
            // Thus we do not need to alter the boundary.
        }
        
        // Compute sigma_{q}
        // (a,b)sigma(k) =
        //   a          if k = sigma^{-1}(b)
        //   b          if k = sigma^{-1}(a)
        //   sigma(k)   else
        // this is done by swapping the values of sigma^{-1}(a) and sigma^{-1}(b) under sigma:
        std::swap( sigma[ sigma_inv[a] ], sigma[ sigma_inv[b] ] );
        
        // sigma^{-1}(a,b) (k) = 
        //   sigma^{-1}(b)  if k = a
        //   sigma^{-1}(a)  if k = b
        //   sigma^{-1}(k)  else
        // this is done by interchanging the values of a and b under sigma^{-1}
        std::swap( sigma_inv[a], sigma_inv[b] );
    }

    // boundary is monotone iff its renormalization is monotone. Thus we check for monotony now to avoid unnecessary renormalization.
    if (not boundary.monotone())
    {
        return Tuple();
    }
    
    // Renormalize all tau'
    for(uint8_t q = 1; q <= boundary.norm(); ++q)
    {
        // Write tau_q = (a,b)
        auto& a = boundary[q].first;
        auto& b = boundary[q].second;

        if(a > k)
        {
            a--;
        }
        if(b > k)
        {
            b--;
        }
    }
    
    boundary.p -= 1;
    
    return boundary;
}

Tuple Tuple :: d_hor_reduced(const uint8_t i) const
{
    const Tuple boundary = this->d_hor(i);
    if( boundary.is_multiple_of_a() == false )
    {
        return boundary;
    }
    else
    {
        return Tuple();
    }
}

Tuple Tuple :: d_hor_double_complex( const uint8_t k ) const
{
    Tuple boundary = *this;
    
    // start with sigma_0.
    Permutation sigma = long_cycle();
    Permutation sigma_inv = long_cycle_inv();
    
    for(uint8_t q = 1; q <= boundary.norm(); ++q)
    {
        // Write tau_q = (a,b)
        auto a = boundary[q].first;
        auto b = boundary[q].second;
        
        //Write (k, sigma_{q-1}(k)) = (k,l)
        auto l = sigma[k];
        
         // Compute tau':
        // Most of the time the transpositions are disjoint hence (a,b)(k,l) = (k,l)(a,b) and
        // the left transposition will be killed by D_k
        if( k != a && k != b && l != a && l != b )
        {
        }
        // The degenerate case:
        // k and l are both part of tau_q.
        else if( a == std::max(k, l) and ( b == std::min(k, l) or k == l))
        {
            return Tuple();
        }
        // The non degenerate case:
        else
        {   
            // Compute Z = (a,b)(k,l)
            // Z(k) = l iff l != a and l != b hence k == a or k == b
            // In this case: (Z(k),k)Z = (k,l)(a,b)(k,l) = (c,l) with c != k
            if( l != a && l != b)
            {
                if(a != k)
                {
                    boundary[q].first  = std::max(a,l);
                    boundary[q].second = std::min(a,l);
                }
                else
                {
                    boundary[q].first  = std::max(b,l);
                    boundary[q].second = std::min(b,l);
                }
            }
            // Z(k) != l iff l = a or l = b hence k != a and k != b
            // In this case: (Z(k),k)Z = (c,k)(a,b)(k,l) = (c,l) with c != l,
            // but we see that this is just (a,b):
            // if l != a, then k,l != a and we map
            // a to a to b != c since otherwise a would map to k
            // therefore a = c and l = b.
            // Thus we do not need to alter the boundary.
        }
        
        // Compute sigma_{q}
        // (a,b)sigma(k) =
        //   a          if k = sigma^{-1}(b)
        //   b          if k = sigma^{-1}(a)
        //   sigma(k)   else
        // this is done by swapping the values of sigma^{-1}(a) and sigma^{-1}(b) under sigma:
        std::swap( sigma[ sigma_inv[a] ], sigma[ sigma_inv[b] ] );
        
        // sigma^{-1}(a,b) (k) = 
        //   sigma^{-1}(b)  if k = a
        //   sigma^{-1}(a)  if k = b
        //   sigma^{-1}(k)  else
        // this is done by interchanging the values of a and b under sigma^{-1}
        std::swap( sigma_inv[a], sigma_inv[b] );
    }
    
    // Renormalize all tau'
    for(uint8_t q = 1; q <= boundary.norm(); ++q)
    {
        // Write tau_q = (a,b)
        auto& a = boundary[q].first;
        auto& b = boundary[q].second;

        if(a > k)
        {
            a--;
        }
        if(b > k)
        {
            b--;
        }
    }
    
    boundary.p -= 1;
    
    return boundary;
}

Permutation Tuple::sigma_h() const
{
    // initialize with sigma_0
    Permutation sigma_inv = long_cycle_inv();

    for (uint8_t i = 1; i <= norm(); ++i)
    {
        // write tau_i = (a, b)
        const uint8_t& a = at(i).first;
        const uint8_t& b = at(i).second;
        
        // compute sigma_i^{-1}
        std::swap( sigma_inv[ a ], sigma_inv[ b ] );
    }
    
    // compute sigma
    Permutation sigma(p+1);
    for( uint8_t i = 0; i <= p; ++i )
    {
        sigma[sigma_inv[i]] = i;
    }
    
    return sigma;
}

std::map< uint8_t, int8_t > Tuple::orientation_sign( ) const
{
    Permutation sigma = sigma_h();
    std::map< uint8_t, Permutation > cycle_decomp = sigma.cycle_decomposition();    
    std::map< uint8_t, int8_t > sign;
    
    uint8_t i = 1; // counter of cycles
    for ( auto it_1 = cycle_decomp.begin(); it_1 != cycle_decomp.end(); ++it_1 )
    {
        Permutation cycle = it_1->second;
        uint8_t min_symbol = it_1->first;
        // if the cycle is a fixpoint (a), we set sign(a) = 0 for the sake of completeness.
        if (cycle.is_fix_point(min_symbol))
        {
            sign[min_symbol] = 0;
            ++i;
            continue;
        }
        uint8_t second_min_symbol = 0;
        // determine the second min symbol of the cycle
        for (size_t m = 1; m < cycle.size(); ++m)
        {
            if (cycle.is_contained(m) && m != min_symbol)
            {
                second_min_symbol = m;
                break;
            }
        }
        // Find k.
        // note that
        //   a_{1,1} < ... < a_{i,1} < b ,
        // hence
        //   a_{i-l,1} < b    for    l >= 0
        // is impossible and we can start to search at the position k = i.
        uint8_t k = i;
        
        // note that since we exclude the case that second_min_sybols will be sorted in at the end,
        // min_symbol can be sorted in between the cycles and this loop will never reach the liniel.
        
        // note that initially (for k = i) it_2.first = a_{k,1} < b and we found our position k iff the first time b < a_{k+1,1}
        // in the other case we have again a_{k+1,1} < b. Induction.
        for ( auto it_2 = (std::next(it_1)); it_2 != cycle_decomp.end(); ++it_2 )
        {
            uint8_t next_min_symbol = it_2->first;
            if ( second_min_symbol < next_min_symbol )
            {
                if ( ((k - i) % 2) == 0 )
                {
                    sign[min_symbol] = 1;
                }
                else
                {
                    sign[min_symbol] = -1;
                }
                break;
            }
            ++k;
        }
        // If b > a_{m, 1}, we still need to determine the sign of min_symbol.
        if (k == cycle_decomp.size())
        {
            if ( ((k - i) % 2) == 0 )
            {
                sign[min_symbol] = 1;
            }
            else
            {
                sign[min_symbol] = -1;
            }
        }
        
        // for all other symbols of the cycle, we set the sign to 1
        for (size_t c = 0; c < cycle.size(); ++c)
        {
            if (cycle.is_contained(c) and c != min_symbol)
            {
                sign[c] = 1;
            }
        }
        ++i;
    }
    return sign;
}

Permutation Tuple::long_cycle() const
{
    Permutation sigma(p+1);
    for(uint8_t k = 0; k < p; ++k)
    {
        sigma[k] = k+1;
    }
    sigma[p] = 0;
    return sigma;
}

Permutation Tuple::long_cycle_inv() const
{
    Permutation sigma(p+1);
    for(uint8_t k = 1; k <= p; ++k)
    {
        sigma[k] = k-1;
    }
    sigma[0] = p;
    return sigma;
}

const std::vector< Transposition >& Tuple::get_data_rep() const
{
    return rep;
}

Tuple create_cell( const size_t h, ... )
{
    Tuple t(h);
    t.p = 0;
    va_list args;
    va_start(args, h);
    
    for ( size_t i = h; i > 0; --i )
    {
        const uint8_t a = va_arg(args, int);
        const uint8_t b = va_arg(args, int);
        t[i] = Transposition( a, b );
        t.p = std::max( t.p, (uint32_t)std::max( a, b ) );
    }
    va_end(args);
    return t;
}

size_t HashTuple :: operator ()( const Tuple &tuple ) const
{
    size_t hashvalue = 0;
    size_t offset = 2;
    for( const auto& cit : tuple.rep )
    {
        hashvalue += offset*(cit.first + 8*cit.second);
        offset *= 16;
    }
    return hashvalue;
}

std::ostream& operator<< (std::ostream& stream, const Tuple& tuple)
{
    for( auto it = tuple.rep.crbegin(); it != tuple.rep.crend(); ++it )
    {
        stream << "(" << (uint32_t)it->first << ",";
        if( it->first != it->second )
        {
            stream << (uint32_t)it->second << ")";
        }
        else
        {
            stream << "*)";
        }
    }
    
    return stream;
}

Tuple operator*( const Tuple& v_2, const Tuple& v_1 )
{
    const uint32_t & p_1 = v_1.p;
    const uint32_t & p_2 = v_2.p;
    
    const size_t h_1 = v_1.norm();
    const size_t h_2 = v_2.norm();
    const size_t h_prod = h_1 + h_2;
    Tuple prod;
    prod.p = p_1 + p_2;
    prod.id = std::numeric_limits< size_t >::max();
    
    auto& rep = prod.rep;
    rep.reserve( h_prod * sizeof( Transposition ) );
    
    for( size_t i = 0; i < h_1; ++i )
    {
        rep.emplace_back( v_1.rep.at(i) );
    }
    for( size_t i = 0; i < h_2; ++i )
    {
        rep.emplace_back( Transposition( p_1 + v_2.rep.at(i).first, p_1 + v_2.rep.at(i).second ) );
    }
    
    return prod;
}
