/* =============================================================================

	A finite volume code to solve the Saint-Venant system in an open channel.

   =============================================================================

	Copyright (C) 2013 Mathieu Besson

	This program is free software: you can redistribute
	it and/or modify it under the terms of the
	GNU General Public License as published by the Free Software Foundation, either 
	version 3 of the License, or (at your option) any later version.

	This program is distributed in the hope that it will be useful, 
	but WITHOUT ANY WARRANTY; without even the implied warranty 
	of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
	See the GNU General Public License for more details.

	You should have received a copy of the GNU General Public License 
	along with this program. If not, see <http://www.gnu.org/licenses/>.*/

/* =============================================================================

		 Bodies of several functions and procedures needed for the code

   ============================================================================= */
#include "usefull.h"
void save(float tab1[nbCell+2],float tab2[nbCell+2],int t)
{
	FILE *fichier1;
	FILE *fichier2;
	int i;
	char nomFichier1[50];
	char nomFichier2[50];

	/* first file for the height. */
	sprintf(nomFichier1,"%d",t);
	strcat(nomFichier1,"_");
	strcat(nomFichier1,"height");
	strcat(nomFichier1,".txt");
	
	/* second file for the discharge. */
	sprintf(nomFichier2,"%d",t);
	strcat(nomFichier2,"_");
	strcat(nomFichier2,"discharge");
	strcat(nomFichier2,".txt");
	
	fichier1 = fopen(nomFichier1,"w");
	for(i=0;i<=nbCell+1;i++)
	{
		fprintf(fichier1,"%f ",tab1[i]);
	}
	fclose(fichier1);
	
	fichier2 = fopen(nomFichier2,"w");
	for(i=0;i<=nbCell+1;i++)
	{
		fprintf(fichier2,"%f ",tab2[i]);
	}
	fclose(fichier2);
}

void saveExact(float tab1[nbCell+2],float tab2[nbCell+2],int t)
{
	FILE *fichier1;
	FILE *fichier2;
	int i;
	char nomFichier1[50];
	char nomFichier2[50];

	/* first file for the height. */
	sprintf(nomFichier1,"%d",t);
	strcat(nomFichier1,"_");
	strcat(nomFichier1,"hExact");
	strcat(nomFichier1,".txt");
	
	/* second file for the discharge. */
	sprintf(nomFichier2,"%d",t);
	strcat(nomFichier2,"_");
	strcat(nomFichier2,"huExact");
	strcat(nomFichier2,".txt");
	
	fichier1 = fopen(nomFichier1,"w");
	for(i=0;i<=nbCell+1;i++)
	{
		fprintf(fichier1,"%f ",tab1[i]);
	}
	fclose(fichier1);
	
	fichier2 = fopen(nomFichier2,"w");
	for(i=0;i<=nbCell+1;i++)
	{
		fprintf(fichier2,"%f ",tab2[i]);
	}
	fclose(fichier2);
}

float maxVal(float u[nbCell+2])
{
	float maxU;
	int i;
	
	maxU = 0;
	for(i=0;i<=nbCell+1;i++)
	{
		if(u[i]>maxU)
		{
			maxU = u[i];
		}
	}
	return maxU;
}	

float max(float a,float b)
{
	if(a>b)
	{
		return a;
	}
	else
	{
		return b;
	}
}

float min(float a,float b)
{
	if(a<b)
	{
		return a;
	}
	else
	{
		return b;
	}
}

float vAbs(float x)
{
	
	float valAbs;
	
	if(x<0)
	{
		valAbs = -x;
	}
	else
	{
		valAbs = x;
	}
	
	return valAbs;
}

void getEigenValues(float hu[nbCell+2],float h[nbCell+2],float EGPlus[nbCell+2],float EGMinus[nbCell+2])
{
	int i;
	
	for(i=0;i<=nbCell+1;i++)
	{
		if(h[i]<=epsilon)
		{
			EGPlus[i] = 0;
			EGMinus[i] = 0;
		}
		else
		{
			EGPlus[i] = hu[i]/h[i] + sqrt(g*h[i]);
			EGMinus[i] = hu[i]/h[i] - sqrt(g*h[i]);
		}
	}
}

float power3Over2(float a) /* a must be non-negative */
{
	float res;
	
	res = a*sqrt(a);
	
	return res;
}

float getNorme(float h[nbCell+2],float hex[nbCell+2],float dx)
{
	float normh = 0;
	float diff[nbCell+2];
	int i;
	
	for(i=0;i<=nbCell+1;i++)
	{
		diff[i] = 0.0;
	}
	
	for(i=1;i<nbCell+1;i++)
	{
		diff[i] = vAbs(h[i] - hex[i]);
	}
	
	normh = maxVal(diff); /* max norm: no dx factor */
	return normh;
}
