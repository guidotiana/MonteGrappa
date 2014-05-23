/*
 * peptides.c
 *
 *  Created on: Apr 6, 2011
 *      Author: guido
 */

#include "montegrappa.h"
#include "grappino.h"

void MakePeptide(struct s_polymer *p, struct s_mc_parms *parms, int from, int to, int wp, int *npol, int npep, struct s_potential *u)
{
	int i,j,k,ip,q=0,o=1;
	double dx,dy,dz;

	for (ip=0;ip<npep;ip++)
	{
		(p+ip+1)->nback = to - from + 1;
		sprintf((p+ip+1)->title,"peptide #%d",ip+1);

		// copy segment of polymer wp as new polymers
		q=0;
		for (i=from;i<=to;i++)
		{
			(((p+ip+1)->back)+q)->ia = (((p+wp)->back)+i)->ia;
			(((p+ip+1)->back)+q)->itype = (((p+wp)->back)+i)->itype;
			strcpy((((p+ip+1)->back)+q)->aa,(((p+wp)->back)+i)->aa);
			strcpy((((p+ip+1)->back)+q)->type,(((p+wp)->back)+i)->type);
			(((p+ip+1)->back)+q)->iaa = (((p+wp)->back)+i)->iaa - from + 1;				// renumber iaa
			((((p+ip+1)->back)+q)->pos).x = ((((p+wp)->back)+i)->pos).x;
			((((p+ip+1)->back)+q)->pos).y = ((((p+wp)->back)+i)->pos).y;
			((((p+ip+1)->back)+q)->pos).z = ((((p+wp)->back)+i)->pos).z;
			((((p+ip+1)->back)+q)->sph).ang = ((((p+wp)->back)+i)->sph).ang;
			((((p+ip+1)->back)+q)->sph).dih = ((((p+wp)->back)+i)->sph).dih;
			((((p+ip+1)->back)+q)->sph).r = ((((p+wp)->back)+i)->sph).r;
			(((p+ip+1)->back)+q)->nside = (((p+wp)->back)+i)->nside;
			(((p+ip+1)->back)+q)->move = (((p+wp)->back)+i)->move;
			(((p+ip+1)->back)+q)->irot = (((p+wp)->back)+i)->irot;
			(((p+ip+1)->back)+q)->nrot = (((p+wp)->back)+i)->nrot;
			(((p+ip+1)->back)+q)->d2_next = (((p+wp)->back)+i)->d2_next;

			for (j=0;j<(((p+wp)->back)+i)->nside;j++)
			{
				strcpy( (((((p+ip+1)->back)+q)->side)+j)->type, (((((p+wp)->back)+i)->side)+j)->type );
				(((((p+ip+1)->back)+q)->side)+j)->itype = (((((p+wp)->back)+i)->side)+j)->itype;
				(((((p+ip+1)->back)+q)->side)+j)->ia = (((((p+wp)->back)+i)->side)+j)->ia;
				((((((p+ip+1)->back)+q)->side)+j)->pos).x = ((((((p+wp)->back)+i)->side)+j)->pos).x;
				((((((p+ip+1)->back)+q)->side)+j)->pos).y = ((((((p+wp)->back)+i)->side)+j)->pos).y;
				((((((p+ip+1)->back)+q)->side)+j)->pos).z = ((((((p+wp)->back)+i)->side)+j)->pos).z;
				for (k=0;k<(((p+wp)->back)+i)->nrot;k++)
				{
					(((((((p+ip+1)->back)+q)->side)+j)->rot)+k)->b1 = (((((((p+wp)->back)+i)->side)+j)->rot)+k)->b1;
					(((((((p+ip+1)->back)+q)->side)+j)->rot)+k)->b2 = (((((((p+wp)->back)+i)->side)+j)->rot)+k)->b2;
					(((((((p+ip+1)->back)+q)->side)+j)->rot)+k)->b3 = (((((((p+wp)->back)+i)->side)+j)->rot)+k)->b3;
					((((((((p+ip+1)->back)+q)->side)+j)->rot)+k)->ang).ang = ((((((((p+wp)->back)+i)->side)+j)->rot)+k)->ang).ang;
					((((((((p+ip+1)->back)+q)->side)+j)->rot)+k)->ang).dih = ((((((((p+wp)->back)+i)->side)+j)->rot)+k)->ang).dih;
					((((((((p+ip+1)->back)+q)->side)+j)->rot)+k)->ang).r = ((((((((p+wp)->back)+i)->side)+j)->rot)+k)->ang).r;
				}
			}

		q ++;
		}
	}

	// shift peptides
	for (i=(*npol);i<(*npol)+npep;i++)
	{
		do
		{
			dx = frand()-0.5; dy = frand()-0.5; dz = frand()-0.5;
			DisplaceCoM(p,i,dx,dy,dz);
			if (EnergyBoxPolymer(p,u,i)) DisplaceCoM(p,i,-dx,-dy,-dz);
			else o = CheckOverlaps(p,u,parms,(*npol)+npep,parms->nosidechains,0);
		} while (o==1);
	}
	(*npol) += npep;

}

