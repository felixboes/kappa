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
