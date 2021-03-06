// The software kappa is a collection of programs to compute the homology of
// the moduli space of surfaces using the radial model.
// Copyright (C) 2013 - 2018  Felix Boes and Anna Hermann
// 
// This file is part of kappa.
// 
// kappa is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// kappa is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with kappa.  If not, see <http://www.gnu.org/licenses/>.


#include "misc.hpp"

bool file_exists( std::string path )
{
    boost::filesystem::path p(path);
    try
    {
        if( boost::filesystem::exists(p) && boost::filesystem::is_regular_file(p) )
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    catch( const boost::filesystem::filesystem_error& ex )
    {
        std::cout << ex.what() << std::endl;
    }
    return false;
}

bool touch( std::string path )
{
    boost::filesystem::path p(path);     
    
    if( file_exists(path) == true )
    {
        try
        {
            boost::filesystem::last_write_time( p , std::time(0) );
            return true;
        }
        catch( const boost::filesystem::filesystem_error& ex )
        {
            std::cout << ex.what() << std::endl;
        }
        return false;
    }
    else
    {
        std::ofstream ofs;
        std::string native_path = p.native();  
        ofs.open( native_path );
        if( ofs.is_open() == false )
        {
            ofs.close();
            return false;
        }
        ofs.close();
        return true;
    }
}

bool directory_exists( std::string path )
{
    boost::filesystem::path p(path);
    try
    {
        if( boost::filesystem::exists(p) && boost::filesystem::is_directory(p) )
        {
           return true;
        }
        else
        {
            return false;
        }
    }    
    catch( const boost::filesystem::filesystem_error& ex )
    {
        std::cout << ex.what() << std::endl;
    }
    return false;
}

bool create_directory( std::string path )
{
    boost::filesystem::path p(path);
    try
    {
        if( boost::filesystem::create_directory(p) == true )
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    catch( const boost::filesystem::filesystem_error& ex )
    {
        std::cout << ex.what() << std::endl;
    }
    return false;
}


bool create_working_directories( bool print_status_messages )
{
    bool status = true;
    std::list<std::string> dirs = {
        "./cache/",
        "./cache/bases_parallel/",
        "./cache/bases_radial/",
        "./cache/differentials_parallel",
        "./cache/differentials_radial",
        "./cache/list_of_files_that_should_not_be_overwritten",
        "./results"
    };
    
    for( const auto& dir : dirs )
    {
        if( directory_exists( dir ) == false )
        {
            if( print_status_messages == true )
            {
                std::cout << "Creating directory '" << dir << "'.";
            }
            if( create_directory( dir ) == false )
            {
                if( print_status_messages == true )
                {
                    std::cout << " Creation failed." << std::endl;
                }
                status = false;
            }
            else
            {
                std::cout << std::endl;
            }
        }
    }
    
    return status;
}

std::string current_date(){
    return boost::posix_time::to_simple_string( boost::posix_time::ptime(boost::posix_time::second_clock::local_time()) );
}

/*
In order to define
  double current_memory_usage_in_mb()
recall the following man page.

man 5 proc

...
proc/[pid]/stat
  Status  information about the process.  This is used by ps(1).  It
  is defined in /usr/src/linux/fs/proc/array.c.

  The fields, in order, with their  proper  scanf(3)  format  speci‐
  fiers, are:

  pid %d      The process ID.

  comm %s     The  filename of the executable, in parentheses.  This
              is visible whether or not the  executable  is  swapped
              out.

  state %c    One character from the string "RSDZTW" where R is run‐
              ning, S is sleeping in an  interruptible  wait,  D  is
              waiting  in uninterruptible disk sleep, Z is zombie, T
              is traced or stopped (on a signal), and W is paging.

  ppid %d     The PID of the parent.

  pgrp %d     The process group ID of the process.

  session %d  The session ID of the process.

  tty_nr %d   The controlling terminal of the process.   (The  minor
              device  number is contained in the combination of bits
              31 to 20 and 7 to 0; the major  device  number  is  in
              bits 15 to 8.)

  tpgid %d    The ID of the foreground process group of the control‐
              ling terminal of the process.

  flags %u (%lu before Linux 2.6.22)
              The kernel flags word of the process.  For  bit  mean‐
              ings,  see the PF_* defines in the Linux kernel source
              file include/linux/sched.h.   Details  depend  on  the
              kernel version.

  minflt %lu  The  number of minor faults the process has made which
              have not required loading a memory page from disk.

  cminflt %lu The number of minor faults that the process's  waited-
              for children have made.

  majflt %lu  The  number of major faults the process has made which
              have required loading a memory page from disk.

  cmajflt %lu The number of major faults that the process's  waited-
              for children have made.

  utime %lu   Amount of time that this process has been scheduled in
              user  mode,  measured  in  clock  ticks   (divide   by
              sysconf(_SC_CLK_TCK)).    This  includes  guest  time,
              guest_time (time spent  running  a  virtual  CPU,  see
              below), so that applications that are not aware of the
              guest time field do not lose that time from their cal‐
              culations.

  stime %lu   Amount of time that this process has been scheduled in
              kernel  mode,  measured  in  clock  ticks  (divide  by
              sysconf(_SC_CLK_TCK)).

  cutime %ld  Amount of time that this process's waited-for children
              have been scheduled in user mode,  measured  in  clock
              ticks  (divide  by  sysconf(_SC_CLK_TCK)).   (See also
              times(2).)   This  includes  guest  time,  cguest_time
              (time spent running a virtual CPU, see below).

  cstime %ld  Amount of time that this process's waited-for children
              have been scheduled in kernel mode, measured in  clock
              ticks (divide by sysconf(_SC_CLK_TCK)).

  priority %ld
              (Explanation  for  Linux  2.6) For processes running a
              real-time  scheduling  policy   (policy   below;   see
              sched_setscheduler(2)), this is the negated scheduling
              priority, minus one; that is, a number in the range -2
              to  -100,  corresponding  to real-time priorities 1 to
              99.   For  processes  running  under  a  non-real-time
              scheduling policy, this is the raw nice value (setpri‐
              ority(2)) as represented in the  kernel.   The  kernel
              stores nice values as numbers in the range 0 (high) to
              39 (low), corresponding to the user-visible nice range
              of -20 to 19.

              Before Linux 2.6, this was a scaled value based on the
              scheduler weighting given to this process.

  nice %ld    The nice value (see setpriority(2)), a  value  in  the
              range 19 (low priority) to -20 (high priority).

  num_threads %ld
              Number  of  threads in this process (since Linux 2.6).
              Before kernel 2.6, this field was hard coded to 0 as a
              placeholder for an earlier removed field.

  itrealvalue %ld
              The time in jiffies before the next SIGALRM is sent to
              the process due to an interval  timer.   Since  kernel
              2.6.17,  this  field  is  no longer maintained, and is
              hard coded as 0.

  starttime %llu (was %lu before Linux 2.6)
              The time the process started after  system  boot.   In
              kernels  before Linux 2.6, this value was expressed in
              jiffies.  Since Linux 2.6, the value is  expressed  in
              clock ticks (divide by sysconf(_SC_CLK_TCK)).

  vsize %lu   Virtual memory size in bytes.
...

*/
double current_memory_usage_in_mb()
{
    int pid;
    std::string comm;
    char state;
    int ppid, pgrp, session, tty_nr, tpgid;
    unsigned int flags;
    unsigned long int minflt, cminflt, majflt, cmajflt, utime, stime;
    long int cutime, cstime, priority, nice, num_threads, itrealvalue;
    long long unsigned int starttime;
    long unsigned int vsize;
    // We do not need more status information.
    
    std::ifstream proc_stat("/proc/self/stat", std::ios_base::in);
    proc_stat >> pid
             >> comm
             >> state
             >> ppid >> pgrp >> session >> tty_nr >> tpgid
             >> flags
             >> minflt >> cminflt >> majflt >> cmajflt >> utime >> stime
             >> cutime >> cstime >> priority >> nice >> num_threads >> itrealvalue
             >> starttime
             >> vsize;
   proc_stat.close();
   return vsize / 1024.0 / 1024.0;
}

void limit_memory( int32_t percent )
{
    long int pages = sysconf(_SC_PHYS_PAGES);
    long int page_size = sysconf(_SC_PAGE_SIZE);
    
    #ifdef __USE_LARGEFILE64
    rlimit64 lim;
    lim.rlim_cur = ( ( ((rlim64_t)pages) * page_size ) / 100 ) * percent;
    lim.rlim_max = RLIM64_INFINITY;
    if( setrlimit64( RLIMIT_AS, &lim ) != 0 )
    {
        std::cout << "An error occured in limit_memory(). " << strerror(errno) << std::endl;
    }
    #else
    rlimit lim;
    lim.rlim_cur = ( ( ((rlim_t)pages) * page_size ) / 100 ) * percent;
    lim.rlim_max = RLIM_INFINITY;
    if( setrlimit( RLIMIT_AS, &lim ) != 0 )
    {
        std::cout << "An error occured in limit_memory(). " << strerror(errno) << std::endl;
    }
    #endif
}

std::string tex_preamble()
{
    std::stringstream tex;
    tex << "\\documentclass[paper=a4, fontsize=11pt, english]{scrreprt}" << std::endl
        << "\\usepackage{amsmath,amssymb,amsthm,amsfonts,amsbsy,latexsym}" << std::endl
        << "\\usepackage{tikz}" << std::endl
        << "\\usepackage[landscape, left=1cm, right=1cm, top=1cm, bottom=1cm]{geometry}" << std::endl
        << "\\linespread{5}" << std::endl
        << "\\setlength{\\parindent}{0pt}" << std::endl
        << "\\begin{document}" << std::endl;
    
    return tex.str();
}

std::string tex_cell( const Tuple& cell )
{
    const int32_t h = cell.norm();
    const int32_t p = cell.p;
    
    std::stringstream tex;
    tex << "\\tikz[baseline={([yshift=-2.5pt]current bounding box.center)}, x=15pt, y=7pt, every node/.style={shape=circle, fill=black, inner sep=.8pt}]{" << std::endl
        << "    \\foreach \\y in {1,...," << p <<"}" << std::endl
        << "    {" << std::endl
        << "        \\draw (-0.5, \\y) -- (" << h-1 << ".5, \\y);" << std::endl
        << "    }" << std::endl
        << "    \\draw[color=black!50] (-0.5,.7) -- (" << h-1 << ".5,.7) -- ("<< h-1 << ".5, "<< p << ".3) -- (-0.5, " << p << ".3) -- (-0.5, .7);" << std::endl;
    for( int32_t i = 1; i <= h; ++i )
    {
    tex << "    \\draw (" << h-i << "," << (int32_t)cell.at(i).first << ") node {} -- (" << h-i << "," << (int32_t)cell.at(i).second << ") node {};" << std::endl;
    }
    tex << "}" << std::endl;
    
    return tex.str();
}

std::string tex_cell( const std::list<Tuple>& cells )
{
    std::stringstream tex;
    for( const auto& it : cells )
    {
        tex << tex_cell( it );
    }
    
    return tex.str();
}

std::string tex_end()
{
    std::stringstream tex;
    tex << "\\end{document}" << std::endl;
    
    return tex.str();
}

std::string kappa_version( int argc, char** argv )
{
    std::stringstream ret;
    if( argc > 0 && argv != nullptr )
    {
        ret << "Programcall:     ";
        for( int i = 0; i < argc; ++i )
        {
            ret << argv[i] << " ";
        }
        ret << std::endl;
    }
    ret << "Program version: " << program_version_by_git << std::endl
        << "GMP version:     " << gmp_version << std::endl
        << "Boost version:   " << BOOST_VERSION / 100000 << "." << BOOST_VERSION / 100 % 1000 << "." << BOOST_VERSION % 100 << std::endl
        << "Date:            " << current_date() << std::endl;
    return ret.str();
}

template<>
std::string filename_prefix_differentials<Q>( const uint32_t g, const uint32_t m, const bool radial )
{
    return std::string("./cache/differentials_") + std::string( radial == true ? "radial" : "parallel" ) + "/q_" + std::to_string(g) + "_" + std::to_string(m) + "_";
}

template<>
std::string filename_prefix_differentials<Zm>( const uint32_t g, const uint32_t m, const bool radial )
{
    return std::string("./cache/differentials_") + std::string( radial == true ? "radial" : "parallel" ) + "/s" + std::to_string( Zm::get_modulus() ) + "_" + std::to_string(g) + "_" + std::to_string(m) + "_";
}

std::vector< std::vector< std::vector< size_t > > > list_set_partitions(size_t n)
{
    std::vector< std::vector< std::vector< size_t > > > res(n);
    
    std::vector< size_t > k( n, 0 );
    std::vector< size_t > m( n, 0 );
    
    auto next_partition = [&]() -> bool
    {
        for( size_t i = n-1; i > 0; --i )
        {
            if( k.at(i) <= m.at(i-1) )
            {
                k[i] += 1;
                if( k.at(i) > m.at(i) )
                {
                    m[i] = k.at(i);
                }
                for( size_t j = i+1; j < n; ++j )
                {
                    k[j] = k[0];
                    m[j] = m[i];
                }
                return true;
            }
        }
        return false;
    };
    

    do
    {
        res[ m.at(n-1) ].push_back( k );
    } while( next_partition() == true );
    
    return res;
}

std::vector< std::vector< std::vector< size_t > > > list_connected_partitions(size_t n)
{
    std::vector< std::vector< std::vector< size_t > > > res(n);
    
    std::vector< size_t > k( n, 0 );
    
    auto next_partition = [&]() -> bool
    {
        for( size_t i = n; i --> 0; )
        {
            k[i] += 1;
            if( k.at(i) < n )
            {
                for( size_t j = i+1; j < n; ++j )
                {
                    k[j] = k[i];
                }
                return true;
            }
        }
        return false;
    };
    
    do
    {
        size_t num_part = 0;
        size_t cur = k.at(0);
        for( size_t i = 1; i < n; ++i )
        {
            if( cur < k.at(i) )
            {
                ++num_part;
                cur = k.at(i);
            }
        }
        res[ num_part ].push_back( k );
    } while( next_partition() == true );
    
    return res;
}
