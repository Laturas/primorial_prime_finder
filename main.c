#include <gmp.h>
#include "primelist.c"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define NEW(i) mpz_t i; mpz_init(i);
#define NEWUI(i,b) mpz_t i; mpz_init_set_ui(i,b)

mpz_t expn;
mpz_t bas;


void primorial_handle(int start, int end, mpz_t primorial, int realbase, mpz_t prime) {
    mpz_set_ui(expn,primes[start]); //printf("expn = %d ", primes[start]);
    for (int k = start + 1; k <= end; k++) {
        mpz_mul_ui(expn,expn,primes[k]);
        //printf("* %d ",primes[k]);
    }
    mpz_divexact(bas,primorial,expn);
    //printf("\nbas = all other than that\n");
    NEWUI(rebase,realbase);
    mpz_powm(bas,rebase,bas,prime);
    mpz_clear(rebase);
}
//#define PRIMECOUNT 3
#define MAX_CHUNKS 32

void lucas_check(const int primorial_index) {
    clock_t start = clock();
    NEWUI(primorial,2);
    for (int i = 1; i <= primorial_index; i++) {
        mpz_mul_ui(primorial,primorial,primes[i]);
    }
    NEW(prime);
    mpz_add_ui(prime,primorial,1);
    gmp_printf("Testing %d#+1 = %Zd\n", primes[primorial_index], prime);
    
    clock_t powbench = 0.0;
    clock_t divbench = 0.0;
    int steps[MAX_CHUNKS] = {0};
    for (int i = 0; i < MAX_CHUNKS; i++) {
        steps[i] = ((primorial_index * (i + 1)) / MAX_CHUNKS) + 1;
    }

    NEW(mpow); NEW(non_one); mpz_init(bas); mpz_init(expn);
    for (int i = 10; i < 1799; i++) {
        mpz_set_ui(mpow,primes[i]);
        mpz_powm(mpow,mpow,primorial,prime);
        if (mpz_cmp_ui(mpow,1) != 0) {printf("a=%d failed initial mod1 test. This is composite\n",i); fflush(stdout); return;}

        mpz_set_ui(mpow,primes[i]);
        int flag = 0;
        int max_steps = MAX_CHUNKS - 2;
        int current_step = 0;
        for (int j = 0; j <= primorial_index; j++) {
            if (j == 0) {primorial_handle(0,steps[0] - 1,primorial,primes[i],prime);}
            if (current_step != max_steps) {
                if (j == steps[current_step]) {primorial_handle(steps[current_step],steps[current_step + 1] - 1,primorial,primes[i],prime); current_step++;}
            }
            else {
                if (j == steps[current_step]) {primorial_handle(steps[current_step],primorial_index,primorial,primes[i],prime);}
            }
            //if (j == steps[0]) {primorial_handle(steps[0],steps[1] - 1,primorial,primes[i],prime);}
            //if (j == steps[1]) {primorial_handle(steps[1],steps[2] - 1,primorial,primes[i],prime);}
            
            mpz_divexact_ui(expn,expn,primes[j]);
            //printf("Removing %d from expn, ",primes[j]);
            if (j != 0 && j != steps[current_step]) {mpz_mul_ui(expn,expn,primes[j-1]); /*printf("adding %d\n",primes[j - 1]);*/}
            //clock_t strt_divbench = clock();
            //mpz_divexact_ui(non_one,primorial,primes[j]);
            //divbench += clock() - strt_divbench;

            //gmp_printf("Testing %Zd ^ %Zd mod %Zd\n", mpow, non_one, prime);
            //clock_t strt_powbench = clock();
            mpz_powm(non_one,bas,expn,prime); // store, base, expn, mod
            //powbench += clock() - strt_powbench;

            flag = (mpz_cmp_ui(non_one,1) == 0);
            if (flag) {printf("a=%d failed on factor %d\n",primes[i],primes[j]); fflush(stdout); break;}
            printf("%d Passed\n",primes[j]);
        }
        if (!flag) {printf("This is prime! Valid found at a = %d\n", primes[i]); break;} // if this broke, I changed the gotos to breaks
        if (i > 100) {printf("BREAKING: Probably composite\n"); break;}
    }
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

void pock_segment(int start, int end, mpz_t bigf, int realbase, mpz_t prime) {
    mpz_set_ui(expn,primes[start]); //printf("expn = %d ", primes[start]);
    for (int k = start + 1; k <= end; k++) {
        mpz_mul_ui(expn,expn,primes[k]);
        //printf("* %d ",primes[k]);
    }
    mpz_divexact(bas,bigf,expn);
    //printf("\nbas = all other than that\n");
    NEWUI(rebase,realbase);
    mpz_powm(bas,rebase,bas,prime);
    mpz_clear(rebase);
}
#undef MAX_CHUNKS
#define MAX_CHUNKS 16

void pocklington_check(const int primorial_index) {
    clock_t start = clock();
    NEWUI(primorial,2);
    for (int i = 1; i <= primorial_index; i++) {
        mpz_mul_ui(primorial,primorial,primes[i]);
    }
    //NEW(prime);
    //mpz_add_ui(prime,primorial,1);
    //gmp_printf("Testing %d#+1 = %Zd\n", primes[primorial_index], prime);
    
    clock_t powbench = 0.0;
    clock_t divbench = 0.0;

    NEW(sqrt); NEWUI(bigf,1); NEW(prime);
    mpz_set(prime,primorial);
    mpz_add_ui(prime,prime,1);
    mpz_sqrt(sqrt, prime);
    //gmp_printf("Square root: %Zd\n",sqrt);
    //gmp_printf("Original: %Zd\n",primorial);
    int bigf_fac_index = -1;
    do {
        bigf_fac_index++;
        mpz_mul_ui(bigf,bigf,primes[primorial_index - bigf_fac_index]);
        //gmp_printf("%Zd\n",bigf);
    } while (mpz_cmp(sqrt,bigf) > 0);
    //printf("Passed after: %d factors\n", j + 1);

    int steps[MAX_CHUNKS] = {0};
    for (int i = 0; i < MAX_CHUNKS; i++) {
        steps[i] = (((primorial_index - bigf_fac_index) * (i + 1)) / MAX_CHUNKS) + 1;
    }

    NEW(mpow); NEW(non_one); mpz_init(bas); mpz_init(expn);
    for (int i = 2; i < 50; i++) {
        mpz_set_ui(mpow,i);
        mpz_powm(mpow,mpow,primorial,prime);
        int current_step = bigf_fac_index;
        int flag = 0;
        if (mpz_cmp_ui(mpow,1) != 0) {printf("a=%d failed initial mod1 test. This is composite\n",i); fflush(stdout); return;}

        const int max_steps = MAX_CHUNKS - 2;
        int item = bigf_fac_index;
        for (int j = 0; item <= primorial_index; j++) {
            item = j + bigf_fac_index;
            if (j == 0) {primorial_handle(bigf_fac_index,bigf_fac_index + steps[0] - 1,primorial,primes[item],prime);}
            if (current_step != max_steps) {
                if (j == steps[current_step]) {primorial_handle(steps[current_step],steps[current_step + 1] - 1,primorial,primes[item],prime); current_step++;}
            }
            else {
                if (j == steps[current_step]) {primorial_handle(steps[current_step],primorial_index,primorial,primes[item],prime);}
            }
            //if (j == steps[0]) {primorial_handle(steps[0],steps[1] - 1,primorial,primes[i],prime);}
            //if (j == steps[1]) {primorial_handle(steps[1],steps[2] - 1,primorial,primes[i],prime);}
            
            mpz_divexact_ui(expn,expn,primes[item]);
            //printf("Removing %d from expn, ",primes[j]);
            if (j != 0 && j != steps[current_step]) {mpz_mul_ui(expn,expn,primes[j-1]); /*printf("adding %d\n",primes[j - 1]);*/}
            //clock_t strt_divbench = clock();
            //mpz_divexact_ui(non_one,primorial,primes[j]);
            //divbench += clock() - strt_divbench;

            gmp_printf("Testing %Zd ^ %Zd mod %Zd\n", bas, expn, prime);
            //clock_t strt_powbench = clock();
            mpz_powm(non_one,bas,expn,prime); // store, base, expn, mod
            //powbench += clock() - strt_powbench;

            flag = (mpz_cmp_ui(non_one,1) == 0);
            if (flag) {printf("a=%d failed on factor %d\n",primes[i],primes[item]); fflush(stdout); break;}
            printf("%d Passed\n",primes[item]); fflush(stdout);
        }
        if (!flag) {printf("This is prime! Valid found at a = %d\n", primes[i]); break;}
        if (i > 100) {printf("BREAKING: Probably composite\n"); break;}
    }

    
    clock_t end = clock();
    double elapsed_time = (end-start)/(double)CLOCKS_PER_SEC;
    printf("Process completed in %lf seconds \n",elapsed_time);
    //printf("powers took: %lf seconds \n",powbench/(double)CLOCKS_PER_SEC);
    //printf("divisions took: %lf seconds \n",divbench/(double)CLOCKS_PER_SEC);
    mpz_clear(sqrt);
    mpz_clear(primorial);
    mpz_clear(bigf);
    //mpz_clear(prime);
}

void chunk_base(int start, int end, mpz_t primorial, int realbase, mpz_t prime) {

}

void ppbls_test(const int primorial_index) {
    clock_t start = clock();
    NEWUI(primorial,2);
    for (int i = 1; i <= primorial_index; i++) {
        mpz_mul_ui(primorial,primorial,primes[i]);
    }
    NEW(sqrt); NEWUI(F,1); NEW(prime);
    mpz_set(prime,primorial);
    mpz_add_ui(prime,prime,1);
    mpz_sqrt(sqrt, prime);
    //gmp_printf("Square root: %Zd\n",sqrt);
    //gmp_printf("Original: %Zd\n",primorial);
    int bigf_fac_index = -1;
    do {
        bigf_fac_index++;
        mpz_mul_ui(F,F,primes[primorial_index - bigf_fac_index]);
        //gmp_printf("%Zd\n",bigf);
    } while (mpz_cmp(sqrt,F) > 0);

    gmp_printf("Testing %d#+1 = %Zd\n", primes[primorial_index], prime);

    NEW(remexp);
    mpz_divexact(remexp,primorial,F);
    int steps[MAX_CHUNKS] = {0};
    for (int i = 0; i < MAX_CHUNKS; i++) {
        steps[i] = (((primorial_index - bigf_fac_index) * (i + 1)) / MAX_CHUNKS) + 1;
    }
    
    NEW(res);
    for (int aui = 2; aui < 3; aui++) {
        printf("Testing witness a = %d^(%d#/F)\n",aui,primes[primorial_index]);
        NEWUI(a,aui);
        mpz_powm(a,a,remexp,prime);
        mpz_powm(res,a,F,prime); // store, base, expn, mod
        int a_flag = (mpz_cmp_ui(res,1) == 0);
        //if (flag) {printf("We win! :D\na = %d\n",aui); fflush(stdout); return;}
        if (!a_flag) {printf("failed for g=%d\n",aui); continue;}
        int current_step = 0;
        int max_steps = MAX_CHUNKS - 2; 
        for (int q = 0; q <= bigf_fac_index; q++) {
            NEW(temppow); NEW(gcd);
            mpz_divexact_ui(temppow,F,primes[primorial_index - q]);
            if (q == 0); //{primorial_handle(bigf_fac_index,bigf_fac_index + steps[0] - 1,primorial,primes[item],prime);}
            if (current_step != max_steps) {
                if (q == steps[current_step]); //{primorial_handle(steps[current_step],steps[current_step + 1] - 1,primorial,primes[item],prime); current_step++;}
            }
            else {
                if (q == steps[current_step]); //{primorial_handle(steps[current_step],primorial_index,primorial,primes[item],prime);}
            }
            //divbench += clock() - strt_divbench;

            //gmp_printf("Testing %Zd ^ %Zd mod %Zd\n", mpow, non_one, prime);
            //clock_t strt_powbench = clock();
            mpz_powm(res,a,temppow,prime); // store, base, expn, mod
            mpz_sub_ui(res,res,1);
            mpz_gcd(gcd,res,prime);
            if ((mpz_cmp_ui(gcd,1) != 0)) {
                mpz_clear(temppow);
                mpz_clear(gcd);
                mpz_clear(a);
                printf("Failed at q = %d\n",primes[primorial_index - q]);
                break;
            }
            printf("q = %d Passed\n",primes[primorial_index - q]);
            mpz_clear(temppow);
            mpz_clear(gcd);
            if (q == bigf_fac_index) {
                printf("%d#+1 is prime! :D\n", primes[primorial_index]);
                //printf("This is prime! :D\n");
                clock_t end = clock();
                double elapsed_time = (end-start)/(double)CLOCKS_PER_SEC;
                printf("Process completed in %lf seconds \n",elapsed_time);
                return;
            }
        }
        mpz_clear(a);
    }
    printf("Probably not prime :(\n");

}

int main(int argc, char* argv[]) {
    int primorial_index = atoi(argv[1]) - 1;
    //lucas_check(primorial_index);
    ppbls_test(primorial_index);
}
