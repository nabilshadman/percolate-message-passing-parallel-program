/* Main header file for percolation code. */

/* System size L */
#define L 768

/* Number of processes NPROC */
#define NPROC 16

/* Prototypes for supplied functions */

/* Visualisation */
void mapwrite(char *percfile, int map[L][L], int ncluster);
void mapwritedynamic(char *percfile, int **map, int l, int ncluster);

/* Random numbers */
void rinit(int ijkl);
float uni(void);
