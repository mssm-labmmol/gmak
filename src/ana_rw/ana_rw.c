#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <omp.h>

#define CHUNK_SIZE 5096
#define MALLOC_VAR(f,n,t) f = (t*) malloc (n * sizeof(t))

#ifdef QUIET
#undef DEBUG
#endif

unsigned int count_lines (char* fn)
{
    unsigned int nlines = 0;
    int bytes;
    char buffer[CHUNK_SIZE + 1], lastchar = '\n';
    FILE* fp = fopen(fn, "r");

    while ((bytes = fread(buffer, sizeof(char), CHUNK_SIZE, fp)))
    {
        if (buffer[0] == '@' || buffer[0] == '#')
            nlines--;
        lastchar = buffer[bytes - 1];
        for (char *c = buffer; (c = memchr(c, '\n', bytes - (c - buffer))); c++)
        {
            if  (c[1] == '@' || c[1] == '#')
                continue;
            nlines++;
        }
    }
    if (lastchar != '\n')
        nlines++;
    fclose(fp);

#ifdef DEBUG
    fprintf(stderr, "debug: file %s has %d lines \n", fn, nlines);
#endif

    return nlines;
}

int read_matrix_from_file (char* fn, int ncols, ...)
{
    char buffer[CHUNK_SIZE + 1];
    int bytes;
    char* line_ptr;
    int data_read = 0;
    /* create array of pointers */
    va_list p_args;
    double** cp = (double**) malloc (ncols*sizeof(double*));
    int offset;
    va_start (p_args, ncols);
    for (int i = 0; i < ncols; ++i)
    {
        cp[i] = va_arg (p_args, double*);
    }
    va_end (p_args);
    FILE* fp = fopen(fn,"r");

    while ( (bytes = fread(buffer, sizeof(char), CHUNK_SIZE, fp) ) )
    /* read chunk */
    {
        /* go back to previous newline character */
        buffer[bytes] = '\0';
        char* point;
        int old_bytes = bytes;
        line_ptr = buffer;
        point = strrchr(buffer, '\n');
        *(point) = '\0';
        fseek(fp, -(buffer+bytes - point - 1)*sizeof(char), SEEK_CUR);
        // TODO this might be wrong
        bytes = ((int) point) - ((int) buffer);
#ifdef DEBUG
            fprintf(stderr, "debug: new chunk:\n[%s]\n", buffer);
            fprintf(stderr, "debug: ignored bytes = %d\n", old_bytes - bytes);
#endif
        /* split into lines */
        do
        {
            int firstchar;
            if (line_ptr >= buffer + bytes - 1)
                break;

            if (line_ptr == buffer)
            {
                firstchar = 0;
            }
            else
            {
                firstchar = 1;
            }

#ifdef DEBUG
            fprintf(stderr, "debug: first char in line is [%c]\n", line_ptr[firstchar]);
#endif
            /* process line */
            if (line_ptr[firstchar] == '@' || line_ptr[firstchar] == '#')
            {
#ifdef DEBUG
                fprintf(stderr, "skipped line in file %s ---\n", fn);
#endif
                line_ptr++;
                continue;
            }
            for (int i = 0; i < ncols; i++)
            {
                if (cp[i] != NULL)
                {
                    sscanf(line_ptr, "%lf%n", &(cp[i][data_read]), &offset);
#ifdef DEBUG
                    fprintf(stderr, "reading col %d for data %d in file %s ---> %lf\n", i, data_read, fn, cp[i][data_read]);
#endif
                } else {
                    double tmp;
                    sscanf(line_ptr, "%lf%n", &tmp, &offset);
#ifdef DEBUG
                    fprintf(stderr, "ignoring col %d for data %d in file %s ---> %lf\n", i, data_read, fn, tmp);
#endif
                }
                line_ptr += offset;
            }
            data_read++;
            //line_ptr++;
        } while ( (line_ptr = memchr(line_ptr, '\n', bytes - (line_ptr - buffer))) );
    }
    fclose(fp);
    return data_read;
}

int boolopt (       int argc,
                 char** argv,
            const char* flag)
{
    for (int i = 0; i < argc; ++i)
    {
        if (!strcmp(flag,argv[i]))
        {
            return 1;
        }
    }
    return 0;
}

char* multiopt2str (int argc,
                 char** argv,
            const char* flag,
            int        entry)
{
    char *out;
    for (int i = 0; i < argc; ++i)
    {
        if (!strcmp(flag,argv[i]))
        {
            int n = strlen(argv[i+entry+1]);
            out = (char*) malloc ((n+1)*sizeof(char));
            strcpy(out, argv[i+1+entry]);
            return out;
        }
    }
    return NULL;
}

double multiopt2double (int argc,
                 char** argv,
            const char* flag,
            int        entry)
{
    char *str;
    double db;
    str = multiopt2str(argc,argv,flag,entry);
    db = atof(str);
    free(str);
    return db;
}

double multiopt2int (int argc,
                 char** argv,
            const char* flag,
            int        entry)
{
    char *str;
    int db;
    str = multiopt2str(argc,argv,flag,entry);
    db = atoi(str);
    free(str);
    return db;
}

int main (int  argc,
        char** argv)
{
    char *fn_tot, *fn_minus_B, *fn_minus_B_sB,
         *fn_minus_B_sB_A, *fn_minus_B_sB_A_sA, *fn_pars;
    char *output_fmt;
    int  ref_idx, npars;
    double *tot, *minus_B, *minus_B_sB, *minus_B_sB_A, *minus_B_sB_A_sA, *pars_A, *pars_B;
    double *baseline, *sA, *A, *sB, *B, *times;
    unsigned int size_tot, size_minus_B, size_minus_B_sB, size_minus_B_sB_A, size_minus_B_sB_A_sA;
    unsigned int num_cores = omp_get_num_procs();
    unsigned int opt_time;

    if (argc == 1)
    {
        fprintf(stdout, "usage: %s -f <0.dat> ... <4.dat> -r <refstate_index> -p <parameters.dat> -o <C_output_format> -time <bool_print_time>\n\n", argv[0]);
        fprintf(stdout, "\t0.dat = total potential energy of trajectory in reference state\n"
                "\t1.dat = above (0.dat) subtracting c12 contributions\n"
                "\t2.dat = above (1.dat) subtracting sqrtc12 constributions\n"
                "\t3.dat = above (2.dat) subtracting c6 contributions\n"
                "\t4.dat = above (3.dat) subtracting sqrtc6 contributions\n");
        return 0;
    }

#ifndef QUIET
    fprintf(stderr, "Using %d OpenMP threads.\n", num_cores);
#endif

    fn_tot = multiopt2str(argc, argv, "-f", 0);
    fn_minus_B = multiopt2str(argc, argv, "-f", 1);
    fn_minus_B_sB = multiopt2str(argc, argv, "-f", 2);
    fn_minus_B_sB_A = multiopt2str(argc, argv, "-f", 3);
    fn_minus_B_sB_A_sA = multiopt2str(argc, argv, "-f", 4);
    ref_idx = multiopt2int(argc,argv,"-r",0);
    fn_pars = multiopt2str(argc, argv, "-p", 0);
    output_fmt = multiopt2str(argc,argv,"-o",0);
    opt_time = boolopt (argc, argv, "-time");

    size_tot = count_lines(fn_tot);
    size_minus_B = count_lines(fn_minus_B);
    size_minus_B_sB = count_lines(fn_minus_B_sB);
    size_minus_B_sB_A = count_lines(fn_minus_B_sB_A);
    size_minus_B_sB_A_sA = count_lines(fn_minus_B_sB_A_sA);
    npars = count_lines(fn_pars);

    if ( size_tot != size_minus_B || size_minus_B != size_minus_B_sB || size_minus_B_sB != size_minus_B_sB_A 
            || size_minus_B_sB_A != size_minus_B_sB_A_sA ) 
    {
        fprintf(stderr, "ERROR: number of lines do not match in input files.\n");
        return 1;
    }

    /* allocate input data */
    times = (double*) malloc (size_tot*sizeof(double));
    /* braces are for local scope in array of variables */
    {
        double** input_data[] =  { &tot, &minus_B, &minus_B_sB, &minus_B_sB_A, &minus_B_sB_A_sA, &baseline, &sA, &A, &sB, &B };
        int      input_sizes[] = { size_tot, size_minus_B, size_minus_B_sB, size_minus_B_sB_A, size_minus_B_sB_A_sA, size_tot, size_tot, size_tot, size_tot, size_tot };
        char*    input_fn[]    = { fn_tot, fn_minus_B, fn_minus_B_sB, fn_minus_B_sB_A, fn_minus_B_sB_A_sA  };
#pragma omp parallel for num_threads(num_cores)
        for (int i = 0; i < 10; i++)
        {
            MALLOC_VAR(*(input_data[i]), input_sizes[i], double);
            if (i < 5)
            {
#ifndef QUIET
                fprintf(stderr, "Reading file %s... ", input_fn[i]);
#endif
                if (opt_time)
                    read_matrix_from_file (input_fn[i], 2, times, *(input_data[i]));
                else
                    read_matrix_from_file (input_fn[i], 2, NULL, *(input_data[i]));
#ifndef QUIET
                fprintf(stderr, "Done.\n");
#endif
            }
        }
    }
    pars_A = (double*) malloc (npars*sizeof(double));
    pars_B = (double*) malloc (npars*sizeof(double));
#ifndef QUIET
    fprintf(stderr, "Reading file %s... ", fn_pars);
#endif
    read_matrix_from_file (fn_pars, 2, pars_A, pars_B);
#ifndef QUIET
    fprintf(stderr, "Done.\n");
#endif

    /* Compute reference components. */
#pragma omp parallel for num_threads(num_cores)
    for (int i = 0; i < size_tot; ++i)
    {
        baseline[i] = minus_B_sB_A_sA[i];
        sA[i] = minus_B_sB_A[i] - baseline[i];
        A[i]  = minus_B_sB[i] - minus_B_sB_A[i];
        sB[i] = minus_B[i] - minus_B_sB[i];
        B[i]  = tot[i] - minus_B[i];
    }

    /* Compute new potentials. */
#pragma omp parallel for num_threads(num_cores)
    for (int i = 0; i < npars; ++i)
    {
#ifndef QUIET
        fprintf(stderr, "Computing potential for state %d/%d.\n", i, npars);
#endif
        /* new totals */
        double* this_total;

        this_total = (double*) malloc (size_tot*sizeof(double));
        if (this_total == NULL)
            fprintf(stderr, "ERROR: allocating memory for state %d \n", i);
        else
#ifndef QUIET
            fprintf(stderr, "success in allocating memory for state %d \n", i);
#endif

        for (int j = 0; j < size_tot; j++)
        {
            this_total[j] = 0.0;
            this_total[j] += baseline[j];
            this_total[j] += sqrt(pars_A[i] / pars_A[ref_idx]) * sA[j];
            this_total[j] += (pars_A[i] / pars_A[ref_idx]) * A[j];
            this_total[j] += sqrt(pars_B[i] / pars_B[ref_idx]) * sB[j];
            this_total[j] += (pars_B[i] / pars_B[ref_idx]) * B[j];
        }
        /* write to file */
#ifndef QUIET
        fprintf(stderr, "Writing state %d to file.\n", i);
#endif
        char this_file_name[512];
        sprintf(this_file_name, output_fmt, i);
        FILE *this_fp = fopen(this_file_name, "w");
        for (int j = 0; j < size_tot; j++)
        {
            if (opt_time)
                fprintf(this_fp, "%18.7lf%18.7lf\n", times[j], this_total[j]); 
            else
                fprintf(this_fp, "%18.7lf\n", this_total[j]); 
        }
        free(this_total);
        fclose(this_fp);
#ifndef QUIET
        fprintf(stderr, "Done state %d.\n", i);
#endif
    }

    /* Free allocated variables. */
    free(fn_tot);
    free(fn_minus_B);
    free(fn_minus_B_sB);
    free(fn_minus_B_sB_A);
    free(fn_minus_B_sB_A_sA);
    free(tot);
    free(minus_B);
    free(minus_B_sB);
    free(minus_B_sB_A);
    free(minus_B_sB_A_sA);
    free(times);
    free(output_fmt);

    return 0;
}
