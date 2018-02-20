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


#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

void signalHandler(int sig) {
    printf("Overflow detected.\n");
    exit(1);
}

int main() {
    signal(SIGILL, &signalHandler);  // for clang++
    signal(SIGABRT, &signalHandler); // for g++ (but ftrapv seems not to work properly)

    int largeInt = INT_MAX;
    int normalInt = 42;
    int overflowInt = largeInt + normalInt;  /* should cause overflow */

    /* if compiling with -ftrapv, we shouldn't get here */
    printf("No overflow detected.\n");
    return 0;
}
