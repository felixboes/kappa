#include "alt_grp_tuple.hpp"

bool AltGrpTuple::radial = false;
uint8_t AltGrpTuple::min_symbol = 1;
uint8_t AltGrpTuple::min_boundary_offset = 1;
uint8_t AltGrpTuple::max_boundary_offset = 1;

Norm2Permutation initialNorm2Permutation(Transposition(0, 0), Transposition(0, 0));

AltGrpTuple::AltGrpTuple() :
    p(0),
    rep()
{
}

AltGrpTuple::AltGrpTuple(const size_t num_entries) :
    p(0),
    rep(num_entries, initialNorm2Permutation)
{
}

AltGrpTuple::AltGrpTuple(const uint8_t symbols, const size_t num_entries) :
    p(symbols),
    rep(num_entries, initialNorm2Permutation)
{
}

Norm2Permutation &AltGrpTuple::operator[](const size_t n)
{
    return rep[n-1];
}

Norm2Permutation &AltGrpTuple::at(const size_t n)
{
    return rep.at(n-1);
}

const Norm2Permutation &AltGrpTuple::at(const size_t n) const
{
    return rep.at(n-1);
}

size_t AltGrpTuple::num_entries() const
{
    return rep.size();
}

bool AltGrpTuple::has_correct_norm() const
{
    for(size_t i=1; i<= num_entries(); i++)
    {
        if( !PermutationManager::has_norm_two( at(i) ) )
        {
            return false;
        }
    }
    return true;
}

bool AltGrpTuple::operator==(const AltGrpTuple &t) const
{
    // Observe how two vectors are compared
    // Operations == and != are performed by first comparing sizes, and if they match, the elements are compared
    // sequentially using algorithm equal, which stops at the first mismatch.
    // Source: http://www.cplusplus.com/reference/vector/vector/operators/
    return ((this->rep == t.rep) and (this->p == t.p));
}

bool AltGrpTuple::operator!=(const AltGrpTuple &t) const
{
    // Observe how two vectors are compared
    // Operations == and != are performed by first comparing sizes, and if they match, the elements are compared sequentially
    // using algorithm equal, which stops at the first mismatch.
    // Source: http://www.cplusplus.com/reference/vector/vector/operators/
    return (this->rep != t.rep);
}

void AltGrpTuple::radial_case()
{
    radial = true;
    min_symbol = 0;
    min_boundary_offset = 0;
    max_boundary_offset = 0;
}

void AltGrpTuple::parallel_case()
{
    radial = false;
    min_symbol = 1;
    min_boundary_offset = 1;
    max_boundary_offset = 1;
}

bool AltGrpTuple::is_radial()
{
    return radial;
}

uint8_t AltGrpTuple::get_min_symbol()
{
    return min_symbol;
}

uint8_t AltGrpTuple::get_min_boundary_offset()
{
    return min_boundary_offset;
}

uint8_t AltGrpTuple::get_max_boundary_offset()
{
    return max_boundary_offset;
}

std::ostream& operator<< (std::ostream& stream, const AltGrpTuple& tuple)
{
    stream << "( ";
    for( auto it = tuple.rep.crbegin(); it != tuple.rep.crend(); ++it )
    {
        stream << "(" << std::to_string(it->first.first) << "," <<  std::to_string(it->first.second) << ")"
               << "(" << std::to_string(it->second.first) << "," <<  std::to_string(it->second.second) << ") ";
    }
    stream << ")";
    return stream;
}

uint8_t AltGrpTuple :: num_cycles() const
{
    // count instead the number of cycles of sigma_out_inv() as this is easier to compute than sigma_out.
    return PermutationManager::num_cycles(sigma_out_inv());
}

bool AltGrpTuple :: has_correct_num_cycles(const size_t m) const
{
    // min_symbol = 1 iff radial = false iff num_components = m + 'exactly one boundary component of the surface'
    return num_cycles() == m + min_symbol;
}

bool AltGrpTuple::write_in_Norm2Notation()
{
    for( size_t i = 1; i <= num_entries(); i++ )
    {
        if( !PermutationManager::write_in_norm2notation(at(i)) )
        {
            return false;
        }
    }
    return true;
}

bool AltGrpTuple::fully_unstable()
{
    if (!write_in_Norm2Notation())
    {
        std::cerr << "Error in 'bool AltGrpTuple::fully_unstable() -> "
                  << "Some Norm2Permutation in the AltGrpTuple does not have norm two." << std::endl;
        return false;
    }
    for( size_t i = 1; i <= num_entries() - 1; i++ )
    {
        if( at(i+1).second.first < at(i).first.first )
        {
            return false;
        }
    }
    return true;
}

bool AltGrpTuple::f(uint8_t i)
{
    if( i == 0 or i >= num_entries())
    {
        std::cerr << "Error in 'bool AltGrpTuple::f( uint8_t i)' -> i="
                  << std::to_string(i) << " and num_entries=" << std::to_string(num_entries())
                  << std::endl;
        return false;
    }

    auto cycle_decomp = PermutationManager::descending_cycle_decomposition_of_product(at(i + 1), at(i));

    //check whether \tau_{i+1}, \tau_i is a geodesic pair
    uint8_t norm_of_product = 0;
    for (size_t j = 0; j < cycle_decomp.size(); ++j)
    {
        norm_of_product += cycle_decomp[j].size() -1;
    }
    if(norm_of_product!=4)
    {
        return false;
    }

    //compute at(i).second
    uint8_t ht = cycle_decomp[0][0];
    uint8_t preimage_of_ht = cycle_decomp[0][cycle_decomp[0].size()-1];
    at(i).second = Transposition(ht, preimage_of_ht);
    cycle_decomp[0].erase(cycle_decomp[0].begin());

    //compute at(i).first
    auto iterator_max_of_zeroth_cycle = std::max_element(cycle_decomp[0].begin(), cycle_decomp[0].end());
    uint8_t max_of_zeroeth_cycle = *iterator_max_of_zeroth_cycle;

    if(cycle_decomp[0].size()>1 and (cycle_decomp.size()==1 or max_of_zeroeth_cycle > cycle_decomp[1][0]))
    {
        ht = max_of_zeroeth_cycle;
        if(cycle_decomp[0][0]==ht)
        {
            preimage_of_ht = cycle_decomp[0][cycle_decomp[0].size()-1];
        }
        else
        {
            preimage_of_ht = *(iterator_max_of_zeroth_cycle -1);
        }
        cycle_decomp[0].erase(iterator_max_of_zeroth_cycle);
    }
    else
    {
        ht = cycle_decomp[1][0];
        preimage_of_ht = cycle_decomp[1][cycle_decomp[1].size()-1];
        cycle_decomp[1].erase(cycle_decomp[1].begin());
    }
    at(i).first = Transposition(ht, preimage_of_ht);

    // erase possible fixed points from cycle_decomp
    if(cycle_decomp.size()>1 and cycle_decomp[1].size()<=1)
    {
        cycle_decomp.erase(cycle_decomp.begin()+1);
    }
    if(cycle_decomp[0].size()<=1)
    {
        cycle_decomp.erase(cycle_decomp.begin());
    }

    //compute at(i+1)
    if(cycle_decomp.size()==2) //then each cycle must be a transposition in at(i+1)
    {
        at(i+1).second = Transposition(cycle_decomp[0][0], cycle_decomp[0][1]);
        at(i+1).first = Transposition(cycle_decomp[1][0], cycle_decomp[1][1]);
        return true;
    }
    else if(cycle_decomp.size()==1) //then the remaining cycle must be a three-cycle
    {
        at(i+1).second = Transposition(cycle_decomp[0][0], cycle_decomp[0][2]);
        at(i+1).first = Transposition(cycle_decomp[0][1], cycle_decomp[0][2]);
        return true;
    }

    std::cerr << "Error in 'bool AltGrpTuple::f(size_t i)' -> Reached impossible case." << std::endl;

    return false;
}

bool AltGrpTuple::phi(uint8_t q, uint8_t i)
{
    if( i == 0 or i > q or q> num_entries())
    {
        std::cerr << "Error in 'bool AltGrpTuple::phi( uint8_t q, uint8_t i)' -> q=" << std::to_string(q)
                  << " and i=" << std::to_string(i) << " and num_entries="
                  << std::to_string(num_entries()) << std::endl;
        return false;
    }

    for( uint8_t j = q-1; j >= i; j-- ) // The loop terminates due to i > 0.
    {
        if(!f(j)) // The norm of the product falls.
        {
            return false;
        }
    }
    write_in_Norm2Notation();
    return true;
}

AltGrpTuple AltGrpTuple::d_hor_double_complex(uint8_t j) const
{
    if(not is_radial() and (j==0 or j >= this->p))
    {
        std::cerr << "Error in 'AltGrpTuple AltGrpTuple::d_hor_double_complex(uint8_t j) const': parallel case with "
                  << "p=" << std::to_string(this->p) << " but j=" << std::to_string(j) << std::endl;
        return AltGrpTuple();
    }
    else if(is_radial() and j > this->p)
    {
        std::cerr << "Error in 'AltGrpTuple AltGrpTuple::d_hor_double_complex(uint8_t j) const': radial case with "
                  << "p=" << std::to_string(this->p) << " but j=" << std::to_string(j) << std::endl;
        return AltGrpTuple();
    }

    AltGrpTuple boundary = *this;
    // initialize sigma_qminus1_at_j for q=1 with sigma_0_at_j
    uint8_t sigma_qminus1_at_j = j+1;
    if(j==this->p)
    {
        sigma_qminus1_at_j = 0;
    }

    for(uint8_t q = 1; q <= boundary.num_entries(); ++q)
    {
        Norm2Permutation tau_q = this->at(q);
        uint8_t sigma_q_at_j = PermutationManager::value_at(tau_q, sigma_qminus1_at_j);
        uint8_t tau_q_at_j = PermutationManager::value_at(tau_q, j);

        if(tau_q_at_j==j)
        {
        }
        else if(sigma_qminus1_at_j==j or sigma_q_at_j==j
                or (sigma_qminus1_at_j!= sigma_q_at_j and tau_q_at_j==sigma_qminus1_at_j) )   //degenerate case
        {
            return AltGrpTuple();
        }
        else if(sigma_qminus1_at_j == sigma_q_at_j)
        {
            PermutationManager::substitute_index(boundary.at(q), j, sigma_q_at_j);
        }
        else if(tau_q_at_j != sigma_qminus1_at_j)
        {
            boundary.at(q).first = Transposition(sigma_q_at_j, sigma_qminus1_at_j);
            boundary.at(q).second = Transposition(sigma_qminus1_at_j, tau_q_at_j);
        }
        else
        {
            std::cerr << "Error in 'AltGrpTuple AltGrpTuple::d_hor_double_complex(uint8_t j) const': "
                  << "reached impossible case" << std::endl;
            return AltGrpTuple();
        }

        sigma_qminus1_at_j = sigma_q_at_j;
    }

    //renormalize boundary by dropping the new fixed point j
    for(uint8_t q = 1; q <= boundary.num_entries(); ++q)
    {
        PermutationManager::drop_fixed_point(boundary.at(q), j);
        PermutationManager::write_in_norm2notation(boundary.at(q));
    }
    boundary.p -= 1;

    return boundary;
}

AltGrpTuple AltGrpTuple::d_hor(uint8_t j) const
{
    AltGrpTuple boundary = d_hor_double_complex(j);
    if (boundary==AltGrpTuple() or not boundary.fully_unstable())
    {
        return AltGrpTuple();
    }
    return boundary;
}

std::map<uint8_t, int8_t> AltGrpTuple::orientation_sign() const
{
    Permutation sigma = sigma_out();
    return PermutationManager::orientation_sign_for_ehrenfried(sigma);
}

Permutation AltGrpTuple::long_cycle() const
{
    return Permutation::long_cycle(p);
}

Permutation AltGrpTuple::long_cycle_inv() const
{
    return Permutation::long_cycle_inv(p);
}

Permutation AltGrpTuple::sigma_out_inv() const
{
    Permutation sigma_inv = long_cycle_inv();

    for (uint8_t i = 1; i <= num_entries(); ++i)
    {
        PermutationManager::multiply(sigma_inv, at(i).second);
        PermutationManager::multiply(sigma_inv, at(i).first);
    }

    return sigma_inv;
}

Permutation AltGrpTuple::sigma_out() const
{
    return PermutationManager::inverse(sigma_out_inv());
}

AltGrpTuple create_alt_grp_tuple( const size_t num_entries, ... )
{
    AltGrpTuple t(num_entries);
    t.p = 0;
    va_list args;
    va_start(args, num_entries);

    for ( size_t i = num_entries; i > 0; --i )
    {
        const uint8_t a = va_arg(args, int);
        const uint8_t b = va_arg(args, int);
        const uint8_t c = va_arg(args, int);
        const uint8_t d = va_arg(args, int);
        t[i] = Norm2Permutation(Transposition(a,b), Transposition(c,d));
        t.p = std::max( t.p, std::max(a,b) );
        t.p = std::max( t.p, std::max(c,d) );
    }
    va_end(args);
    return t;
}

size_t HashAltGrpTuple :: operator ()( const AltGrpTuple &tuple ) const
{
    size_t hashvalue = 0;
    size_t offset = 2;
    for( const auto& cit : tuple.rep )
    {
        hashvalue += offset*(cit.first.first + 4*cit.first.second + 16*cit.second.first + 64*cit.second.second);
        offset *= 256;
    }
    return hashvalue;
}

bool AltGrpTuple::has_no_common_fixed_points_except_zero()
{
    std::vector< bool > contained_in_support(p + 1);
    for(uint8_t perm_index=1; perm_index <= num_entries(); perm_index++)
    {
        Norm2Permutation t = at(perm_index) ;
        contained_in_support.at(t.first.first)=true;
        contained_in_support.at(t.first.second)=true;
        contained_in_support.at(t.second.first)=true;
        contained_in_support.at(t.second.second)=true;
    }

    for(uint8_t symbol = 1; symbol <=p; symbol++)
    {
        if(not contained_in_support[symbol])
        {
            return false;
        }
    }
    return true;
}