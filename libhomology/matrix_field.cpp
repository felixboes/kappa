#include "matrix_field.hpp"

template<>
void save_to_file_bz2( const MatrixField<Q>& matrix, std::string filename, const bool print_duration )
{
    if( print_duration == true )
    {
        std::cout << "Saving '" + filename + ".bz2'";
        std::cout.flush();
    }
    Clock measure_duration;

    // Open bz2 stream.
    std::string bz2_command_string = "bzip2 -c -z -7 > " + filename + ".bz2";
    const char* bz2_command = bz2_command_string.c_str();
    FILE* bz2_pipe;
    
    if( (bz2_pipe = popen(bz2_command, "w")) == nullptr )
    {
        std::cout << "Error: Could not open '" << filename << "'." << std::endl;
        return;
    }
   
    // Store number of rows and number of columns.
    fprintf( bz2_pipe, "%zu %zu\n", matrix.num_rows, matrix.num_cols );
    
    // Store number of diagonal entries as well as all entries.
    fprintf( bz2_pipe, "%zu\n", matrix.diagonal.size() );
    for( const auto& entry : matrix.diagonal )
    {
        fprintf( bz2_pipe, "%zu %zu\n", entry.first, entry.second );
    }
    
    // Store matrix entries.
    for( size_t t = 0; t < matrix.num_rows * matrix.num_cols; ++t )
    {
        const auto& entry = matrix.data[t];
        if( mpz_out_raw( bz2_pipe, entry.get_num_mpz_t() ) == 0 )
        {
            std::cout << "Error: Could write numerator. Closing file." << std::endl;
            pclose(bz2_pipe);
        }
        if ( mpz_out_raw( bz2_pipe, entry.get_den_mpz_t() ) == 0 )
        {
            std::cout << "Error: Could write denominator. Closing file." << std::endl;
            pclose(bz2_pipe);
        }
    }
    
    pclose(bz2_pipe);
    
    if( print_duration == true )
    {
        std::cout << " done. Duration: " << measure_duration.duration() << " seconds." << std::endl;
        std::cout.flush();
    }
}

template <>
MatrixField<Q> load_from_file_bz2( std::string filename, const bool print_duration )
{
    if( print_duration == true )
    {
        std::cout << "Loading '" + filename + ".bz2'";
        std::cout.flush();
    }
    Clock measure_duration;

    // Open bz2 stream.
    std::string bz2_command_string = "bzip2 -c -d < " + filename + ".bz2";
    const char* bz2_command = bz2_command_string.c_str();
    FILE* bz2_pipe;
    
    if( (bz2_pipe = popen(bz2_command, "r")) == nullptr )
    {
        std::cout << "Error: Could not open '" << filename << "'." << std::endl;
        return MatrixField<Q>();
    }
    
    // Load number of rows and columns.
    size_t num_rows = 0;
    size_t num_cols = 0;
    if( fscanf( bz2_pipe, "%i %i\n", &num_rows, &num_cols ) == 0 )
    {
        std::cout << "Error: Could not read number of columns or rows. Closing file." << std::endl;
        return MatrixField<Q>();
    }
    MatrixField<Q> m ( num_rows, num_cols );
    
    // Load number of diagonal entries as well as all entries.
    size_t diagonal_size = 0;
    size_t i = 0;
    size_t j = 0;
    auto& diag = m.diagonal;
    if( fscanf( bz2_pipe, "%zu\n", &diagonal_size ) == 0 )
    {
        std::cout << "Error: Could not read number of diagonal entries. Closing file." << std::endl;
        return MatrixField<Q>();
    }
    for( size_t t = 0; t < diagonal_size; ++t )
    {
        if( fscanf( bz2_pipe, "%zu %zu\n", &i, &j ) == 0 )
        {
            std::cout << "Error: Could not read diagonal entries. Closing file." << std::endl;
            return MatrixField<Q>();
        }
        diag.emplace_back(i,j);
    }
    
    // Load matrix entries.
    for( size_t t = 0; t < num_rows * num_cols; ++t )
    {
        auto& entry = m.data[t];
        if( mpz_inp_raw( entry.get_num_mpz_t(), bz2_pipe ) == 0 )
        {
            std::cout << "Error: Could read numerator. Closing file." << std::endl;
            pclose(bz2_pipe);
            return m;
        }
        if( mpz_inp_raw( entry.get_den_mpz_t(), bz2_pipe ) == 0 )
        {
            std::cout << "Error: Could read numerator. Closing file." << std::endl;
            pclose(bz2_pipe);
            return m;
        }
    }
    
    pclose(bz2_pipe);
    
    if( print_duration == true )
    {
        std::cout << " done. Duration: " << measure_duration.duration() << " seconds." << std::endl;
        std::cout.flush();
    }
      
    return m;
}

MatrixBool::MatrixBool() : diagonal(), data(), num_rows(0), num_cols(0)
{
}

MatrixBool::MatrixBool( size_t number_rows, size_t number_cols ) : diagonal(), data( number_rows, boost::dynamic_bitset<>(number_cols, 0) ), num_rows(number_rows), num_cols(number_cols)
{
}

void MatrixBool::row_operation( size_t row_1, size_t row_2, size_t )
{
    data[row_2] ^= data[row_1];
}

 bool MatrixBool::operator()( size_t i, size_t j )
 {
     return data[i][j];
 }

 void MatrixBool::add_entry( size_t i, size_t j)
 {
        data[i][j].flip();
 }

 void MatrixBool::resize (size_t size1, size_t size2, bool )
 {
     for (size_t i = 0; i < data.size(); ++i)
     {
         data[i].reset();
     }

     data.resize(size1);
     for (size_t i = 0; i < size1; ++i)
     {
         data[i].resize(size2);
     }
     num_rows = size1;
     num_cols = size2;
     
     diagonal.clear();
 }
 
 bool MatrixBool::at(size_t i, size_t j) const
{
    return data[i][j];
}

size_t MatrixBool::size1() const
{
    return num_rows;
}

size_t MatrixBool::size2() const
{
    return num_cols;
}

void MatrixBool::clear()
{
    for( auto & it : data )
    {
        it.clear();
    }
}

void MatrixBool :: print_base_changes_in_short_form() const
{
    if( diagonal.size() == 0 )
    {
        std::cout << "The matrix seems to be not diagonalized: The diagonal of the matrix is empty." << std::endl;
        return;      
    }
    
    for( size_t i = 0; i < num_rows; ++i )
    {
        auto diag_entry = diagonal.cbegin();
        // print elements befor the diagonal entry.
        while( diag_entry != diagonal.cend() && diag_entry->first != i )
        {
            std::cout << std::setw(3) << (int32_t)this->at( i, diag_entry->second );
            ++diag_entry;
        }
        // print diagonal entry.
        if( diag_entry != diagonal.cend() )
        {
            std::cout << std::setw(3) << (int32_t)this->at( diag_entry->first, diag_entry->second );
            ++diag_entry;
        }
        // print elements after diagonal entry.
        while( diag_entry != diagonal.cend() )
        {
            std::cout << std::setw(3) << 0;
            ++diag_entry;
        }
        std::cout << std::endl;
    }
}

void MatrixBool :: print_triangular_shape() const
{
    if( diagonal.size() == 0 )
    {
        std::cout << "The matrix seems to be not diagonalized: The diagonal of the matrix is empty." << std::endl;
        return;      
    }
    
    // prepare fast access to the rows storing a diagonal entry.
    std::vector< bool >   diagonal_entry_occures_in_row (num_rows, false);
    std::vector< size_t > diagonal_entry_col_row (num_rows, 0);
    
    for( const auto& diag_entry : diagonal )
    {
        diagonal_entry_occures_in_row[diag_entry.first] = true;
        diagonal_entry_col_row[diag_entry.first] = diag_entry.second;
    }
    
    // print trinagular shape.
    for( size_t i = 0; i < num_rows; ++i )
    {
        // row with diagonal entry.
        if( diagonal_entry_occures_in_row[i] == true )
        {
            size_t j = 0;
            // zeros befor diagonal entry.
            while( j < num_cols && j < diagonal_entry_col_row[i] )
            {
                std::cout << std::setw(3) << "0";
                ++j;
            }
            while( j < num_cols ) // diagonal entry and elements to the right.
            {
                std::cout << std::setw(3) << (int32_t)this->at(i,j);
                ++j;
            }
        }
        else // row of zeros.
        {
            for( size_t j = 0; j < num_cols; ++j )
            {
                std::cout << std::setw(3) << "0";
            }
        }
        std::cout << std::endl;
    }
}

void MatrixBool :: cache_matrix( std::string ) const
{
    std::cout << "Error: MatrixBool cannot be cached since boost::dynamic_bitset is not supported by boost::serialization." << std::endl
              << "       We refere to https://svn.boost.org/trac/boost/ticket/3328" << std::endl;
    return;
}

void MatrixBool :: cache_base_change( std::string ) const
{
    std::cout << "Error: MatrixBool cannot be cached since boost::dynamic_bitset is not supported by boost::serialization." << std::endl
              << "       We refere to https://svn.boost.org/trac/boost/ticket/3328" << std::endl;
    return;
}

void MatrixBool :: cache_triangular_shape( std::string ) const
{
    std::cout << "Error: MatrixBool cannot be cached since boost::dynamic_bitset is not supported by boost::serialization." << std::endl
              << "       We refere to https://svn.boost.org/trac/boost/ticket/3328" << std::endl;
    return;
}

void MatrixBool :: cache_diagonal( std::string filename ) const
{
    // Observe that
    //     save_to_file_bz2 ( diagonal, filename );
    // wont compile because of 'Argument-dependent name lookup'
    save_to_file_bz2( this->diagonal, filename ) ;
}

std::ostream& operator<< ( std::ostream& stream, const MatrixBool & matrix)
{
    for( size_t i = 0; i < matrix.num_rows; ++i )
    {
        for( size_t j = 0; j < matrix.num_cols; )
        {
            stream << matrix.at(i,j);
            if( ++j < matrix.num_cols )
            {
                stream << ",";
            }
            else
            {
                stream << std::endl;
            }
        }
    }
    return stream;
}

/////////////////////////////////////////////////////////////////////////////////////////////////

MatrixBoolCSS::MatrixBoolCSS() : diagonal(), data(), sec_data(), num_rows(0), num_cols(0), sec_num_rows(0), sec_num_cols(0)
{
}

MatrixBoolCSS::MatrixBoolCSS( const MatrixBoolCSSInitialization ini, const size_t num_rows1, const size_t num_cols1, const size_t num_rows2, const size_t num_cols2 ) :
    diagonal(), 
    data(), sec_data(),
    num_rows(0), num_cols(0), sec_num_rows(0), sec_num_cols(0)
{
    if( ini == only_main )
    {
        data = MatrixStorageType( num_rows1, boost::dynamic_bitset<>(num_cols1, 0) );
        num_rows = num_rows1;
        num_cols = num_cols1;
    }
    else if( ini == only_secondary )
    {
        sec_data = MatrixStorageType( num_rows1, boost::dynamic_bitset<>(num_cols1, 0) );
        sec_num_rows = num_rows1;
        sec_num_cols = num_cols1;
    }
    else if( ini == both )
    {
        data = MatrixStorageType( num_rows1, boost::dynamic_bitset<>(num_cols1, 0) );
        num_rows = num_rows1;
        num_cols = num_cols1;
        sec_data = MatrixStorageType( num_rows2, boost::dynamic_bitset<>(num_cols2, 0) );
        sec_num_rows = num_rows2;
        sec_num_cols = num_cols2;
    }
}

MatrixBoolCSS::MatrixBoolCSS( const size_t number_rows, const size_t number_cols ) :
    data( number_rows, boost::dynamic_bitset<>(number_cols, 0) ),
    sec_data(),
    num_rows(number_rows),
    num_cols(number_cols),
    sec_num_rows(0),
    sec_num_cols(0)
{
}

void MatrixBoolCSS::define_operations( const OperationType rt )
{
    if( rt == main_and_secondary )
    {
        row_operation_funct = &ThisType::row_operation_main_and_secondary;
        op_funct = &ThisType::main_op;
        at_funct = &ThisType::main_at;
        add_entry_funct = &ThisType::main_add_entry;
        size1_funct = &ThisType::main_size1;
        size2_funct = &ThisType::main_size2;
    }
    else
    {
        row_operation_funct = &ThisType::row_operation_secondary;
        op_funct = &ThisType::sec_op;
        at_funct = &ThisType::sec_at;
        add_entry_funct = &ThisType::sec_add_entry;
        size1_funct = &ThisType::sec_size1;
        size2_funct = &ThisType::sec_size2;
    }
}

bool MatrixBoolCSS::operator()( const size_t i, const size_t j )
{
   return (this->*op_funct)( i, j );
}


bool MatrixBoolCSS::main_op( const size_t i, const size_t j )
{
    return data[i][j];
}


bool MatrixBoolCSS::sec_op( const size_t i, const size_t j )
{
    return sec_data[i][j];
}


bool MatrixBoolCSS::at( const size_t i, const size_t j ) const
{
   return (this->*at_funct)( i, j );
}


bool MatrixBoolCSS::main_at(const size_t i, const size_t j) const
{
    return data[i][j];
}


bool MatrixBoolCSS::sec_at(const size_t i, const size_t j) const
{
    return sec_data[i][j];
}

void MatrixBoolCSS::row_operation( const size_t row_1, const size_t row_2, const size_t col )
{
   (this->*row_operation_funct)(row_1, row_2, col);
}

void MatrixBoolCSS::row_operation_main_and_secondary( const size_t row_1, const size_t row_2, const size_t )
{
    data[row_2] ^= data[row_1];
    sec_data[row_2] ^= sec_data[row_1];
}

void MatrixBoolCSS::row_operation_secondary( const size_t row_1, const size_t row_2, const size_t )
{
    sec_data[row_2] ^= sec_data[row_1];
}

void MatrixBoolCSS::add_entry( const size_t i, const size_t j )
{
    (this->*add_entry_funct)(i, j);
}

void MatrixBoolCSS::main_add_entry(const size_t i, const size_t j)
{
    data[i][j].flip();
}

void MatrixBoolCSS::sec_add_entry(const size_t i, const size_t j)
{
    sec_data[i][j].flip();
}

void MatrixBoolCSS::main_set( const size_t i, const size_t j, const bool val )
{
    data[i].set(j,val);
}

void MatrixBoolCSS::sec_set( const size_t i, const size_t j, const bool val )
{
    sec_data[i].set(j,val);
}

MatrixBoolCSS::MatrixRowType& MatrixBoolCSS::main_row( const size_t i )
{
    return data[i];
}

const MatrixBoolCSS::MatrixRowType& MatrixBoolCSS::main_row_at( const size_t i ) const
{
    return data.at(i);
}

MatrixBoolCSS::MatrixRowType& MatrixBoolCSS::sec_row( const size_t i )
{
    return sec_data[i];
}

const MatrixBoolCSS::MatrixRowType& MatrixBoolCSS::sec_row_at( const size_t i ) const
{
    return sec_data.at(i);
}

void MatrixBoolCSS::resize(size_t size1, size_t size2, bool)
{
    for(size_t i = 0; i < data.size(); ++i)
    {
        data[i].reset();
    }
    
    data.resize(size1);
    for (size_t i = 0; i < size1; ++i)
    {
        data[i].resize(size2);
    }
    num_rows = size1;
    num_cols = size2;
}

void MatrixBoolCSS::sec_resize (const size_t size1, const size_t size2, const bool)
{
    for(size_t i = 0; i < sec_data.size(); ++i)
    {
        sec_data[i].reset();
    }
    
    sec_data.resize(size1);
    for (size_t i = 0; i < size1; ++i)
    {
        sec_data[i].resize(size2);
    }
    sec_num_rows = size1;
    sec_num_cols = size2;
}

size_t MatrixBoolCSS::size1() const
{
    return (this->*size1_funct)();
}

size_t MatrixBoolCSS::size2() const
{
    return (this->*size2_funct)();
}

size_t MatrixBoolCSS::main_size1() const
{
    return num_rows;
}

size_t MatrixBoolCSS::main_size2() const
{
    return num_cols;
}

size_t MatrixBoolCSS::sec_size1() const
{
    return sec_num_rows;
}

size_t MatrixBoolCSS::sec_size2() const
{
    return sec_num_cols;
}

void MatrixBoolCSS::clear()
{
    for( auto & it : data )
    {
        it.clear();
    }
}

void MatrixBoolCSS::sec_clear()
{
    for( auto & it : sec_data )
    {
        boost::dynamic_bitset<> & row = it;
        for ( size_t j = 0; j < sec_num_cols; ++j)
        {
            row[j] = 0;
        }
    }
}

std::ostream& operator<< ( std::ostream& stream, const MatrixBoolCSS & matrix )
{
    stream << "main " << matrix.num_rows << "x" << matrix.num_cols << ", secondary " << matrix.sec_num_rows << "x" << matrix.sec_num_cols << std::endl;
    stream << "main:" << std::endl;
    for( size_t i = 0; i < matrix.num_rows; ++i )
    {
        for( size_t j = 0; j < matrix.num_cols; )
        {
            stream << matrix.main_at(i,j);
            if( ++j < matrix.num_cols )
            {
                stream << ",";
            }
            else
            {
                stream << std::endl;
            }
        }
    }
    
    stream << "secondary:" << std::endl;
    for( size_t i = 0; i < matrix.sec_num_rows; ++i )
    {
        for( size_t j = 0; j < matrix.sec_num_cols; )
        {
            stream << matrix.sec_at(i,j);
            if( ++j < matrix.sec_num_cols )
            {
                stream << ",";
            }
            else
            {
                stream << std::endl;
            }
        }
    }
    return stream;
}
