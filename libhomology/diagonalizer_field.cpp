#include "diagonalizer_field.hpp"
#include "diagonalizer_field_impl.ipp"

#include "homology.hpp"

// For bool versions of the matrices, we do not have an efficient base change algorithm yet.
template<>
void DiagonalizerField< MatrixBool >::apply_base_changes( MatrixType& differential, const MatrixType& base_changes )
{
    (void)differential;
    (void)base_changes;
}

template<>
void DiagonalizerField< MatrixBoolCSS >::apply_base_changes( MatrixType& differential, const MatrixType& base_changes )
{
    (void)differential;
    (void)base_changes;
}


#define force_template_instantiation(MatrixType)\
    template class DiagonalizerField<MatrixType>;

force_template_instantiation(MatrixQ)
force_template_instantiation(MatrixZm)
force_template_instantiation(MatrixBool)
force_template_instantiation(MatrixQCSS)
force_template_instantiation(MatrixZmCSS)
force_template_instantiation(MatrixBoolCSS)

#undef force_template_instantiation

