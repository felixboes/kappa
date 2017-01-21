#include "monocomplex.hpp"
#include "monocomplex_impl.ipp"


template<>
void update_differential(MatrixBool &           differential,
                         const size_t           row,
                         const size_t           column,
                         const int32_t          ,
                         const int8_t           ,
                         const int8_t           ,
                         const SignConvention & )
{
    differential.add_entry(row, column);
}

template<>
void update_differential(MatrixBoolCSS &        differential,
                         const size_t           row,
                         const size_t           column,
                         const int32_t          ,
                         const int8_t           ,
                         const int8_t           ,
                         const SignConvention & )
{
    differential.add_entry(row, column);
}

/* Force template instantiation for used types */

// We do not store homchains for Zm right now.
template<>
void MonoComplex< ChainComplexZm > :: homchain(int32_t p, bool homology, int32_t maxdimension)
{
    (void)p;
    (void)homology;
    (void)maxdimension;
}
#define force_template_instantiation(MatrixComplex) \
    template class MonoComplex<MatrixComplex>;\
    template void update_differential(MatrixComplex &differential, const size_t row, const size_t column, const int32_t, const int8_t, const int8_t, const SignConvention &);\
    template void monocomplex_work(MonoComplex<MatrixComplex> & monocomplex, MonocomplexWork & work, const uint32_t p, MatrixComplex::MatrixType & differential);

force_template_instantiation(ChainComplexQ)
force_template_instantiation(ChainComplexZm)
force_template_instantiation(ChainComplexZStorageOnly)

#undef force_template_instantiation\

int32_t sign(const int32_t          parity,
             const int8_t           i,
             const int8_t           or_sign,
             const SignConvention & sign_conv )
{
    if ( sign_conv == no_signs)
    {
        return 1;
    }
    if( sign_conv == all_signs )
    {
        int32_t actual_parity = (parity + i) % 2;
        if ( or_sign == -1 )
        {
            actual_parity = (actual_parity + 1) % 2;
        }
        //std::cout << it << " " << i << ": The d^hor_i boundary of " << current_basis << ". This is " << boundary << std::endl;
        //std::cout << it.id << "->" << boundary.id << " in " << "M_{" << basis_complex[p-1].size() << "," << basis_complex[p].size() << "} parity=" << actual_parity << std::endl;
        //std::cout << std::endl;
        if ( actual_parity == 0 )
        {
            return 1;
        }
        else
        {
            return -1;
        }
    }
    else if( sign_conv == no_orientation_sign )
    {
        if ( (parity + i) % 2 == 0 )
        {
            return 1;
        }
        else
        {
            return -1;
        }
    }

    return 0;
}

