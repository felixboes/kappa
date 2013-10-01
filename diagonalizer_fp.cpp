#include "diagonalizer_fp.hpp"

DiagonalizerFp::DiagonalizerFp(MatrixFp &pre_differential, MatrixFp &post_differential) :
    pre(pre_differential),
    post(post_differential)
{
    // Compute the dimension of the kernel of pre.
    // This done by computing the number of lineary independant columns or rows.

    // Is #rows < #cols?
    if( pre.size1() < pre.size2() )
    {
        // Um Vertauschungen zu vermeiden:
        // Iteriere durch einfach verkettete Liste der Zeilen.
        // Lösche Eintrag, falls es in dieser Zeile (bei gegebener Spalte) ein invertierbares Element gibt.
        // Mache Zeilenop auf alle vorhandenen Einträge.
    }
}
