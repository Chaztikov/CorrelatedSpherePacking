#include <stdio.h>

#define MAX2D 1000
#define MAX3RD 1000
#define MAXPIX 1000

FILE *from, *to;

struct porestr {
	unsigned phase: 1;
	unsigned fstart: 10;
	unsigned fend: 10;
	unsigned connct: 10;
} pore[MAXPIX][MAXPIX];

void get_poredat(int poredat[MAX2D][MAX2D], int dim);
int readdim(void);


void main(void)
{
   int poredat[MAX2D][MAX2D];
	int dim, nz;

	/*open file for input*/

	if((from=fopen("pore.out","r"))==NULL)
	{
      printf("Cannot open <pore.out> file.\n");
		exit(1);
	}

	/*open files for output*/

	if((to=fopen("porestr.out","w"))==NULL)
	{
      printf("Cannot open <porestr.out> file.\n");
		exit(1);
	}

	/*read dimension of pore network*/

	dim = readdim();
	printf("The dimension for pore network is : %d\n",dim);

   /*get features of 2D pore*/

   for(nz=0; nz<dim; nz++)
	{
	   get_poredat(poredat, dim);
				
	}

}

int readdim(void)
{
	char dummystr[20];
   int dim;

   fscanf(from, "%s %s %s %s %s %s %d", dummystr, dummystr,
			 dummystr, dummystr, dummystr, dummystr, &dim);
	return dim;
}

void get_poredat(int poredat[][MAX2D], int dim)
{
	int i, j, temp;

   for(j=0; j<dim; j++)
	{
		for(i=0; i<dim; i++)
		{
         fscanf(from, "%d %d %d %d", &temp, &temp, &temp, &temp);
			poredat[j][i] = temp;
      }
    }
}

int get_feature(int dim, int nz, int poredat[][MAX2D], 
					 struct porestr[][MAXPIX])
{
	int col, count, phase;
	int i, j, maxf_col;

   col = count = maxf_col = 0;

   printf("Enter the phase (0 or 1): ");
	scanf("%d", &phase);

   for(j=0; j<dim; j++)
	{
      for(i=0; i<dim; i++)
		{
			temp = poredat[j][i];
			if(temp==phase)
			{
				count++;
				porestr[nz][col].fstart = (unsigned) i;
			}
			while(

		}
	}

