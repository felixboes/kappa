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
        // Wir gehen davon aus, dass die Anzahl der Zeilen echt kleiner als die Anzahl der Spalten ist.
        //
        // Um Dimension des Kerns und des Bildes zu berechnen, reicht es aus,
        // die Anzahl der linear unabhängigen Zeilen zu bestimmen.
        //
        // Wir können im Stufenform-Algorithmus Vertauschungen zu vermeiden:
        // - Iteriere durch die Spalten.
        //   - Iteriere durch einfach verkettete Liste der Zeilen.
        //   - Lösche Eintrag, falls es in dieser Zeile (bei gegebener Spalte) ein invertierbares Element gibt.
        //   - Mache Zeilenop auf alle vorhandenen Einträge.
        size_t row_c = pre.size2();
        for( size_t j = 0; j < row_c; j++ )
        {

        }
    }
}
