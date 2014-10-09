#include <cstdio>
#include <iostream>
#include <string>

#include <gmpxx.h>

/*
--------------------------
man 3 popen:
--------------------------

Synopsis
--------------------------

#include <stdio.h>

FILE *popen(const char *command, const char *type);

int pclose(FILE *stream);

Feature Test Macro Requirements for glibc (see feature_test_macros(7)):

popen(), pclose():
    _POSIX_C_SOURCE >= 2 || _XOPEN_SOURCE || _BSD_SOURCE || _SVID_SOURCE 


Description
--------------------------

The popen() function opens a process by creating a pipe, forking, and invoking the shell. Since a pipe is by definition unidirectional, the type argument may specify only reading or writing, not both; the resulting stream is correspondingly read-only or write-only.

The command argument is a pointer to a null-terminated string containing a shell command line. This command is passed to /bin/sh using the -c flag; interpretation, if any, is performed by the shell. The type argument is a pointer to a null-terminated string which must contain either the letter 'r' for reading or the letter 'w' for writing. Since glibc 2.9, this argument can additionally include the letter 'e', which causes the close-on-exec flag (FD_CLOEXEC) to be set on the underlying file descriptor; see the description of the O_CLOEXEC flag in open(2) for reasons why this may be useful.

The return value from popen() is a normal standard I/O stream in all respects save that it must be closed with pclose() rather than fclose(3). Writing to such a stream writes to the standard input of the command; the command's standard output is the same as that of the process that called popen(), unless this is altered by the command itself. Conversely, reading from a "popened" stream reads the command's standard output, and the command's standard input is the same as that of the process that called popen().

Note that output popen() streams are fully buffered by default.

The pclose() function waits for the associated process to terminate and returns the exit status of the command as returned by wait4(2).
Return Value

The popen() function returns NULL if the fork(2) or pipe(2) calls fail, or if it cannot allocate memory.

The pclose() function returns -1 if wait4(2) returns an error, or some other error is detected.


Errors
--------------------------

The popen() function does not set errno if memory allocation fails. If the underlying fork(2) or pipe(2) fails, errno is set appropriately. If the type argument is invalid, and this condition is detected, errno is set to EINVAL.

If pclose() cannot obtain the child status, errno is set to ECHILD.

--------------------------
man 1 bzip2
--------------------------
.....

OPTIONS
       -c --stdout
              Compress or decompress to standard output.

       -d --decompress
              Force decompression.  bzip2, bunzip2 and bzcat are really the same program, and the decision about what actions to take is done on the basis of which
              name is used.  This flag overrides that mechanism, and forces bzip2 to decompress.

       -z --compress
              The complement to -d: forces compression, regardless of the invocation name.
...
       -1 (or --fast) to -9 (or --best)
              Set  the  block  size  to  100 k, 200 k ...  900 k when compressing.  Has no effect when decompressing.  See MEMORY MANAGEMENT below.  The --fast and
              --best aliases are primarily for GNU gzip compatibility.  In particular, --fast doesn't make things significantly faster.  And --best merely  selects
              the default behaviour.
...

--------------------------
gmp documentation
--------------------------
... 

- Function: size_t mpz_out_raw (FILE *stream, const mpz_t op)
    Output op on stdio stream stream, in raw binary format. The integer is written in a portable format, with 4 bytes of size information, and that many bytes of limbs. Both the size and the limbs are written in decreasing significance order (i.e., in big-endian).
    The output can be read with mpz_inp_raw.
    Return the number of bytes written, or if an error occurred, return 0.
    The output of this can not be read by mpz_inp_raw from GMP 1, because of changes necessary for compatibility between 32-bit and 64-bit machines. 

- Function: size_t mpz_inp_raw (mpz_t rop, FILE *stream)
    Input from stdio stream stream in the format written by mpz_out_raw, and put the result in rop. Return the number of bytes read, or if an error occurred, return 0.
    This routine can read the output from mpz_out_raw also from GMP 1, in spite of changes necessary for compatibility between 32-bit and 64-bit machines. 

...

*/
int main( int argc, char** argv )
{
    if( argc == 1 || std::string(argv[1]) == std::string("c") )
    {
        std::string bz2_command_string = std::string("bzip2 -c -z -7 > test.bz2");
        const char* bz2_command = bz2_command_string.c_str();
        FILE* bz2_pipe;
        
        if( (bz2_pipe = popen(bz2_command, "w")) == nullptr )
        {
            std::cout << "Fehler beim Öffnen." << std::endl;
            return 1;
        }
        
        fprintf( bz2_pipe, "%i %i\n", 40, 24 );
        
        mpq_class t=200;
        t /= 7;
        
        std::cout << t << std::endl;
        
        mpz_out_raw( bz2_pipe, t.get_num().get_mpz_t() );
        mpz_out_raw( bz2_pipe, t.get_den().get_mpz_t() );
        
        t*= mpq_class("234892734987324/8589439837");
        std::cout << t << std::endl;
        
        mpz_out_raw( bz2_pipe, t.get_num().get_mpz_t() );
        mpz_out_raw( bz2_pipe, t.get_den().get_mpz_t() );
        
        pclose(bz2_pipe);
        return 0;
    }
    else
    {
        std::string bz2_command_string = std::string("bzip2 -c -d < test.bz2");
        const char* bz2_command = bz2_command_string.c_str();
        FILE* bz2_pipe;
        
        if( (bz2_pipe = popen(bz2_command, "r")) == nullptr )
        {
            std::cout << "Fehler beim Öffnen." << std::endl;
            return 1;
        }
        
        mpz_t in;
        int a = 0;
        int b = 0;
        fscanf( bz2_pipe, "%i %i\n", &a, &b );
        std::cout << a << " " << b << std::endl;
        mpz_init(in);
        while( mpz_inp_raw( in, bz2_pipe ) != 0 )
        {
            std::cout << in << std::endl;
        }
        mpz_clear(in);
        
        pclose(bz2_pipe);
        return 0;
    }
}
