#include <gmp.h>
#include "primelist.c"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define NEW(i) mpz_t i; mpz_init(i);
#define NEWUI(i,b) mpz_t i; mpz_init_set_ui(i,b)

//#define PRIMECOUNT 3

int main(int argc, char* argv[]) {
    int primorial_index = atoi(argv[1]) - 1;
    clock_t start = clock();
    NEWUI(primorial,2);
    for (int i = 1; i <= primorial_index; i++) {
        mpz_mul_ui(primorial,primorial,primes[i]);
    }
    mpz_t prime;
    mpz_add_ui(prime,primorial,1);
    gmp_printf("Testing %d#+1 = %Zd\n", primes[primorial_index], prime);
    
    clock_t powbench = 0.0;
    clock_t divbench = 0.0;

    NEW(mpow); NEW(non_one);
    for (int i = 3; i < 1799; i++) {
        mpz_set_ui(mpow,primes[i]);
        mpz_powm(mpow,mpow,primorial,prime);
        if (mpz_cmp_ui(mpow,1) != 0) {printf("a=%d failed initial mod1 test\n",i); fflush(stdout); continue;}

        mpz_set_ui(mpow,primes[i]);
        int flag = 0;
        for (int j = 0; j <= primorial_index; j++) {
            clock_t strt_divbench = clock();
            mpz_divexact_ui(non_one,primorial,primes[j]);
            divbench += clock() - strt_divbench;

            //gmp_printf("Testing %Zd ^ %Zd mod %Zd\n", mpow, non_one, prime);
            clock_t strt_powbench = clock();
            mpz_powm(non_one,mpow,non_one,prime); // store, base, exp, mod
            powbench += clock() - strt_powbench;

            flag = (mpz_cmp_ui(non_one,1) == 0);
            if (flag) {printf("a=%d failed on factor %d\n",primes[i],primes[j]); fflush(stdout); break;}
            printf("%d Passed\n",primes[j]);
        }
        if (!flag) {printf("This is prime! Valid found at a = %d\n", primes[i]); goto EXIT;}
    }
    EXITNPRIME:
    printf("Not prime :(\n");
    EXIT:
    clock_t end = clock();
    double elapsed_time = (end-start)/(double)CLOCKS_PER_SEC;
    printf("Process completed in %lf seconds \n",elapsed_time);
    //printf("powers took: %lf seconds \n",powbench/(double)CLOCKS_PER_SEC);
    //printf("divisions took: %lf seconds \n",divbench/(double)CLOCKS_PER_SEC);
    mpz_clear(non_one);
    mpz_clear(primorial);
    mpz_clear(mpow);
    mpz_clear(prime);
}