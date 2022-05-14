/* Simple parallel program to test for percolation of a cluster */
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "percolate.h"

int main(int argc, char *argv[])
{

    /* Variables that define the simulation */
    int seed;
    double rho;

    /* Local variables */
    int M, N;
    int i, j, nhole, step, maxstep, oldval, newval;
    int nchangelocal, nchange, printfreq;
    int itop, ibot, perc;
    double squares, mapsumlocal, mapsumglobal, mapavg;
    double r;

    /* MPI variables */
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm cartcomm;

    int size, rank;
    int tag = 1;
    int reorder, cartrank;
    double tstart, tstop, ttotal, tperstep;

    int cartcoords[2];
    int neighboursranks[4];

    enum DIRECTIONS
    {
        LEFT,
        RIGHT,
        DOWN,
        UP
    };

    /* Initialise MPI, and compute size and rank */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    /* Abort program if NPROC does not equal size */
    if (NPROC != size)
    {
        if (rank == 0)
        {
            printf("percolate: ERROR, NPROC = %d but running on %d\n",
                   NPROC, size);
        }

        MPI_Finalize();
        return 0;
    }

    /* Abort program if command line arguments are missing */
    if (argc < 2 || argc > 3)
    {
        if (rank == 0)
        {
            printf("Usage: percolate <seed> [rho]\n");
        }

        MPI_Finalize();
        return 0;
    }

    /*
     * Set random number seed.
     * Set most important value: the rock density rho (between 0 and 1).
     */
    if (argc == 2)
    {
        seed = atoi(argv[1]);
        rho = 0.4040;
    }
    else if (argc == 3)
    {
        seed = atoi(argv[1]);
        rho = atof(argv[2]);
    }

    /* Create two-dimensional virtual cartesian grid topology */

    /* Decompose processes in 2D cartesian grid */
    int dims[2] = {0, 0};
    MPI_Dims_create(size, 2, dims);

    /* Make horizontal (i-th) dimension periodic */
    int periods[2] = {1, 0};

    /* Create communicator given 2D grid topology */
    reorder = 0;
    MPI_Cart_create(comm, 2, dims, periods, reorder, &cartcomm);

    /* Compute rank in new communicator */
    MPI_Comm_rank(cartcomm, &cartrank);

    /* Compute coordinates in new communicator */
    MPI_Cart_coords(cartcomm, cartrank, 2, cartcoords);

    /* Compute neighbours */

    /* Consider dims[0] = X */
    MPI_Cart_shift(cartcomm, 0, 1, &neighboursranks[LEFT], &neighboursranks[RIGHT]);

    /* Consider dims[1] = Y */
    MPI_Cart_shift(cartcomm, 1, 1, &neighboursranks[DOWN], &neighboursranks[UP]);

    /*
     * Use 2D decomposition over NPROC processes across both dimensions.
     * For an LxL simulation, the local arrays are of size MxN.
     * Set M and N values based on values returned from MPI_Dims_create().
     * M is the row dimension, and N is the column dimension. Each process
     * is ensured to have the same values of M and N for having load balance.
     *
     * Note that the map needs to divide equally amongst all processes, i.e.
     * the values of M and N need to be equal for all processes. If either
     * M or N does not divide L, an error is returned and the program aborts,
     * indicating the user to run the program with another value of L.
     */
    M = L / dims[0];
    N = L / dims[1];

    /* Abort program if M does not divide L */
    if (L % M != 0)
    {
        if (rank == 0)
        {
            printf("percolate: ERROR, M (i.e. row dimension) does not divide L. Try a different value for L.\n");
        }
        MPI_Finalize();
        return 0;
    }

    /* Abort program if N does not divide L */
    if (L % N != 0)
    {
        if (rank == 0)
        {
            printf("percolate: ERROR, N (i.e. column dimension) does not divide L. Try a different value for L.\n");
        }
        MPI_Finalize();
        return 0;
    }

    /* Define the main arrays for the simulation */
    int old[M + 2][N + 2], new[M + 2][N + 2];

    /*
     *  Additional array WITHOUT halos for initialisation and IO. This
     *  is of size LxL because, even in our parallel program, we do
     *  these two steps in serial.
     */
    int map[L][L];
    int maptmp[L][L];

    /* Array to store local part of map */
    int smallmap[M][N];

    /* Set maxstep */
    maxstep = 5 * L;

    if (rank == 0)
    {
        printf("percolate: running on %d process(es)\n", size);

        printf("percolate: L = %d, rho = %f, seed = %d, maxstep = %d\n",
               L, rho, seed, maxstep);

        /* Initialise generator */
        rinit(seed);

        /*
         *  Initialise map with density rho. Zero indicates rock, a positive
         *  value indicates a hole. For the algorithm to work, all the holes
         *  must be initialised with a unique integer.
         */
        nhole = 0;

        for (i = 0; i < L; i++)
        {
            for (j = 0; j < L; j++)
            {
                r = uni();

                if (r < rho)
                {
                    map[i][j] = 0;
                }
                else
                {
                    nhole++;
                    map[i][j] = nhole;
                }
            }
        }

        printf("percolate: rho = %f, actual density = %f\n",
               rho, 1.0 - ((double)nhole) / ((double)L * L));
    }

    /*
     * Use broadcast and copy-back to distribute the map. This is not as
     * elegant as using scatter in the 1D decomposition, but generalises
     * to a 2D decomposition (while scatter does not). Use &map[0][0]
     * syntax as this also works for dynamically allocated arrays.
     */
    MPI_Bcast(&map[0][0], L * L, MPI_INT, 0, comm);

    /*
     * Copy the appropriate section back to smallmap. Could probably
     * eliminate use of smallmap in its entirety, but leave it in here
     * for simplicity.
     */
    for (i = 0; i < M; i++)
    {
        for (j = 0; j < N; j++)
        {
            smallmap[i][j] = map[cartcoords[0] * M + i][cartcoords[1] * N + j];
        }
    }

    /*
     * Initialise the old array: copy the LxL array smallmap to the centre of
     * old, and set the halo values to zero.
     */
    for (i = 1; i <= M; i++)
    {
        for (j = 1; j <= N; j++)
        {
            old[i][j] = smallmap[i - 1][j - 1];
        }
    }

    /* Zero the bottom and top halos */
    for (i = 0; i <= M + 1; i++)
    {
        old[i][0] = 0;
        old[i][N + 1] = 0;
    }

    /* Zero the left and right halos */
    for (j = 0; j <= N + 1; j++)
    {
        old[0][j] = 0;
        old[M + 1][j] = 0;
    }

    /* Initialise core update loop variables*/
    step = 1;
    nchange = 1;
    mapsumglobal = 0.0;
    squares = ((double)L * L);
    printfreq = 100;

    /* Define variables necessary for non-blocking communications */

    /* Create arrays to store requests and statuses for non-blocking halo swaps */
    MPI_Request requests[4];
    MPI_Status statuses[4];

    /* Create vector to communicate halos up and down */
    MPI_Datatype columntype;
    MPI_Type_vector(M, 1, N + 2, MPI_INT, &columntype);
    MPI_Type_commit(&columntype);

    /* Create arrays to receive halos from up and down */
    int halosup[M];
    int halosdown[M];

    /* Set stride and block for non-periodic region in horizontal dimension */
    int jpbcstride = 10;
    int jpbcblock = 8;
    int jremain;

    /* Start timing core update loop */
    MPI_Barrier(comm);
    tstart = MPI_Wtime();

    /*
     * Start core percolate algorithm. Update grid until max step
     * is reached or number of changes is zero.
     */
    while (step <= maxstep)
    {
        /*
         * Swap halos left, right, up, and down.
         * Communications is done using non-blocking communications (i.e.
         * combination of issend/recv in this case).
         */

        /* Call non-blocking synchronous send routines */
        MPI_Issend(&old[1][1], N, MPI_INT, neighboursranks[LEFT], tag,
                   cartcomm, &requests[LEFT]);

        MPI_Issend(&old[M][1], N, MPI_INT, neighboursranks[RIGHT], tag,
                   cartcomm, &requests[RIGHT]);

        MPI_Issend(&(old[1][1]), 1, columntype, neighboursranks[DOWN], tag,
                   cartcomm, &requests[DOWN]);

        MPI_Issend(&(old[1][N]), 1, columntype, neighboursranks[UP], tag,
                   cartcomm, &requests[UP]);

        /*
         * Call corresponding blocking receives for the non-blocking sends
         * by the neighbouring processes
         */
        MPI_Recv(&old[M + 1][1], N, MPI_INT, neighboursranks[RIGHT], tag,
                 cartcomm, &statuses[RIGHT]);

        MPI_Recv(&old[0][1], N, MPI_INT, neighboursranks[LEFT], tag,
                 cartcomm, &statuses[LEFT]);

        MPI_Recv(halosup, M, MPI_INT, neighboursranks[UP], tag,
                 cartcomm, &statuses[UP]);

        MPI_Recv(halosdown, M, MPI_INT, neighboursranks[DOWN], tag,
                 cartcomm, &statuses[DOWN]);

        /* Wait for all non-blocking sends to complete */
        MPI_Waitall(4, requests, statuses);

        /* Copy halos received down/up to their appropriate locations in old array */
        for (i = 1; i <= M; i++)
        {
            old[i][0] = halosdown[i - 1];
            old[i][N + 1] = halosup[i - 1];
        }

        /*
         *  Map is mostly periodic across the first dimension "i".
         *
         *  However, it is only periodic for certain elements: in
         *  blocks of "jpbcblock" cells separated by "jpbcstride". For
         *  example, if jpcblock=4 and jpbcstride=6 then boundary
         *  elements 1 and 2 will be non-periodic, 3, 4, 5 and 6
         *  periodic, 7 and 8 non-periodic, 9, 10, 11 and 12 periodic,
         *  etc. etc.
         *
         *  The simplest approach is to copy across the whole boundary
         *  then discard the unwanted (non-periodic) elements.  It may
         *  seem wasteful to copy the halos then zero some of them, but
         *  this approach is probably the simplest to implement in
         *  parallel.
         */
        for (j = 1; j <= N; j++)
        {
            /* Check if is this is one of the non-periodic elements */
            jremain = (j - 1) % jpbcstride;

            if (jremain < (jpbcstride - jpbcblock))
            {
                /* Periodic boundaries do not operate in this region */
                old[0][j] = 0;
                old[M + 1][j] = 0;
            }
        }

        nchangelocal = 0;

        for (i = 1; i <= M; i++)
        {
            for (j = 1; j <= N; j++)
            {
                oldval = old[i][j];
                newval = oldval;

                /*
                 * Set new[i][j] to be the maximum value of old[i][j]
                 * and its four nearest neighbours
                 */
                if (oldval != 0)
                {
                    if (old[i][j - 1] > newval)
                        newval = old[i][j - 1];
                    if (old[i][j + 1] > newval)
                        newval = old[i][j + 1];
                    if (old[i - 1][j] > newval)
                        newval = old[i - 1][j];
                    if (old[i + 1][j] > newval)
                        newval = old[i + 1][j];

                    if (newval != oldval)
                    {
                        ++nchangelocal;
                    }
                }

                new[i][j] = newval;
            }
        }

        /* Compute global number of changes */
        MPI_Allreduce(&nchangelocal, &nchange, 1, MPI_INT, MPI_SUM, comm);

        /* Report progress every now and then */
        if (step % printfreq == 0)
        {
            /* Compute average value of map on rank 0 */
            mapsumlocal = 0.0;

            for (i = 1; i <= M; i++)
            {
                for (j = 1; j <= N; j++)
                {
                    mapsumlocal = mapsumlocal + ((double)new[i][j]);
                }
            }

            MPI_Reduce(&mapsumlocal, &mapsumglobal, 1, MPI_DOUBLE, MPI_SUM, 0, comm);

            /* Print number of changes and map average */
            if (rank == 0)
            {
                mapavg = mapsumglobal / squares;
                printf("percolate: step %d - changes = %d\npercolate: step %d - map average = %.3f\n",
                       step, nchange, step, mapavg);
            }
        }

        /* Stop computation if number of changes is zero */
        if (nchange == 0)
        {
            /* Compute average value of map on rank 0 */
            mapsumlocal = 0.0;

            for (i = 1; i <= M; i++)
            {
                for (j = 1; j <= N; j++)
                {
                    mapsumlocal = mapsumlocal + ((double)new[i][j]);
                }
            }

            MPI_Reduce(&mapsumlocal, &mapsumglobal, 1, MPI_DOUBLE, MPI_SUM, 0, comm);

            /* Print number of changes and map average */
            if (rank == 0)
            {
                mapavg = mapsumglobal / squares;
                printf("percolate: step %d - changes = %d\npercolate: step %d - map average = %.3f\n",
                       step, nchange, step, mapavg);
                printf("percolate: algorithm completed\n");
            }

            break;
        }

        /* Copy back in preparation for next step, omitting halos */
        for (i = 1; i <= M; i++)
        {
            for (j = 1; j <= N; j++)
            {
                old[i][j] = new[i][j];
            }
        }

        step++;
    }

    /* Stop timing core update loop */
    MPI_Barrier(comm);
    tstop = MPI_Wtime();

    /* Print timing data */
    ttotal = tstop - tstart;
    tperstep = ttotal / ((double)step);

    if (rank == 0)
    {
        printf("percolate: total time to update grid = %.6f second(s), time per step = %.6f second(s)\n",
               ttotal, tperstep);
    }

    /*
     *  We set a maximum number of steps to ensure the algorithm always
     *  terminates. However, if we hit this limit before the algorithm
     *  has finished then there must have been a problem (e.g. the value
     *  of maxstep is too small)
     */
    if (rank == 0)
    {
        if (nchange != 0)
        {
            printf("percolate: WARNING max steps = %d reached but nchange != 0\n",
                   maxstep);
        }
    }

    /* Copy the centre of old, excluding the halos, into smallmap */
    for (i = 1; i <= M; i++)
    {
        for (j = 1; j <= N; j++)
        {
            smallmap[i - 1][j - 1] = old[i][j];
        }
    }

    /*
     *  Use copy and reduce to collect the map.  This is not as elegant
     *  as using gather in the 1D decomposition, but generalises to a 2D
     *  decomposition (while gather does not).  Use &map[0][0] syntax in
     *  reduce as this also works for dynamically allocated arrays.
     */

    /* Zero maptmp */
    for (i = 0; i < L; i++)
    {
        for (j = 0; j < L; j++)
        {
            maptmp[i][j] = 0;
        }
    }

    /* Copy smallmap to correct place in maptmp */
    for (i = 0; i < M; i++)
    {
        for (j = 0; j < N; j++)
        {
            maptmp[cartcoords[0] * M + i][cartcoords[1] * N + j] = smallmap[i][j];
        }
    }

    MPI_Reduce(&(maptmp[0][0]), &(map[0][0]), L * L, MPI_INT, MPI_SUM, 0, comm);

    /*
     *  Test to see if percolation occurred by looking for positive numbers
     *  that appear on both the top and bottom edges
     */
    if (rank == 0)
    {
        perc = 0;

        for (itop = 0; itop < L; itop++)
        {
            if (map[itop][L - 1] > 0)
            {
                for (ibot = 0; ibot < L; ibot++)
                {
                    if (map[ibot][0] == map[itop][L - 1])
                    {
                        perc = 1;
                    }
                }
            }
        }

        if (perc != 0)
        {
            printf("percolate: cluster DOES percolate\n");
        }
        else
        {
            printf("percolate: cluster DOES NOT percolate\n");
        }

        /*
         *  Write the map to the file "map.pgm", displaying the two
         *  largest clusters. If the last argument here was 3, it would
         *  display the three largest clusters etc. The picture looks
         *  cleanest with only a single cluster, but multiple clusters
         *  are useful for debugging.
         */
        mapwrite("map.pgm", map, 2);
    }

    /* Finalise MPI */
    MPI_Finalize();

    return 0;
}
