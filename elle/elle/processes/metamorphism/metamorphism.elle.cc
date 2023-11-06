#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include "nodes.h"
#include "display.h"
#include "check.h"
#include "errnum.h"
#include "error.h"
#include "general.h"
#include "runopts.h"
#include "file.h"
#include "interface.h"
#include "polygon.h"
#include "stats.h"
#include "init.h"
#include "crossings.h"
#include "convert.h"
#include "lut.h"
#include "timefn.h"
#include "mineraldb.h"
#include "update.h"
#include "metamorphism.h"
 #include <algorithm>

using std::vector;

/********************physical constants*********************************/
double energyofdislocations;        // in Jm-1    (from mineraldb)
double energyofsurface;         // in Jm-2    (from mineraldb)
double mobilityofboundary;        // in m2s-1J-1    (aka fudge factor)
double dislocationdensityscaling=10e13; // in m-2    (ie 1.5 in elle file = 1.5e13 in real world)
double truetimestep = 3.1536e10;    // in s     (1000yrs here, not including leap years)
double lengthscale = 1e-3;        // in m     (sides of box)
double R = 8.314472;    // in Jmol-1K-1    (Gas constant)
double Qgbm = 200;        // in Jmole-1    (Activation energy for GBM  (made up number))
/***********************************************************************/

using std::cout;
using std::endl;

//#define DEBUG


double  GetAngleEnergy(double angle);
int Phase_MoveNode(int n, Coords *movedir, int *same,double areaphase,
                 double factor, double fudge,double areaequil,
                 double energyphasephase, double energyhostphase,
                 double energyhosthost);
int GetBodyNodeEnergy(int n, double *total_energy);
int GetCSLFactor(Coords_3D xyz[3], float *factor);
int GetSurfaceNodeEnergy(int n,int meltmineral,Coords *xy,double *energy,
                         double areaphase,double fudge, double areaequil,
                         double energyphasephase, double energyhostphase,
                         double energyhosthost);
double Calcareaphase(int meltmineral);
void ElleNodeSameAttrib(int node, int *same, int *samemin, int attrib_id);
extern double ElleSwitchLength();

double    incroffset=0.002; /* this is reassigned in GBE_GrainGrowth */
int min_names[6]={0,5,4,1,2,3}; // switch mineral names to match modes.txt order (Grt	Opx	Cpx	Pl	q )
double modal[1000][6]; // data from modes.txt
double modal_out[6]; // linear interpolation of modal data for current temperature
double surface_energy[5][3]={ .05,.05,0.05,
							  .05,.05,0.05,
							  .05,.05,0.05,
							  .05,.05,0.05,
							  .05,.05,0.05};
int no_pts; // no of rows in modes.txt
int Metamorphism();
int InitMet();

UserData CurrData; // data read in from command line
int PhaseId; // Phase to be equilabrated


int InitMet()
{
    int err=0;
    char *infile;
	long tmp;
	FILE *modes;
	int i;

	tmp = currenttime(); // seed for random number generator
	#ifdef DEBUG
	tmp = 1028731772; // uncomment define of DEBUG to allow specified random number seed
	#endif
	cout << "seed= " << tmp << endl;
    srand(tmp);

    ElleReinit();
    ElleSetRunFunction(Metamorphism);

    infile = ElleFile();
    if (strlen(infile)>0)
    {
        if (err=ElleReadData(infile)) OnError(infile,err);


    }
    else /* warn that there is no file open */;

    ElleUserData(CurrData);

	/*
	 * read in modal data for const P series.
	 * It would probably be nicer to generate PT modal grids for each phase so we can wander aimlessly around PT space
	 */
	modes=fopen("modes.txt","r");

	i=0;
	do
	{
		err=fscanf(modes,"%lf %lf %lf %lf %lf %lf",&modal[i][0],&modal[i][1],&modal[i][2],&modal[i][3],&modal[i][4],&modal[i][5]);
		i++;
	}while(err> 0);

	fclose(modes);
	no_pts=i-1;

}
/*!
 \brief	principal calculation of multiple phase reactions based on input thermodynaic data from Perplex
 \par	Description:

 		This calculation is based on the melt code but includes the possibility of evolving multiple phases based on the
		calculated equilibrium modal % based on Perplex calculations, which at present have to be done before hand,
		The code loops through each phase indivisually to try and move two-phase boundaries so as to try and minimise the area
		misfit compared to the equilibrium area. As present the surface energies are set by an internal surface energy array, however
		this should be externalised. The PID control using the fudge factor also needs to be made more robust. The behaviour of this system
		as a phase is about to disappear is not very stable and needs to be modified so that multiple smaller timesteps
		replace the single large one. The method of defining which phases are present is all very ad hoc, but we need to define the
		nature of relationship with Perplex before we can improve this.

		usage elle_metamorphism -i filename -s number of time steps -f interval to save elle files -u mobility_of_grain_boundaries fudge_factor incremental_Temp incremental_Pressure

		e.g.

		elle_metamorphism -i eqm_900.elle -s 300 -f 25 -u 1e-10 .001 0.5 0


 */

int Metamorphism()
{
    int i, j, k, n;			// indices
	int  same;				// same phase across boundaries?
	int max;				// max node number (not max number of nodes)
    FILE *modal_stats;		// output file pointer
    Coords incr; 			// calculated incremental node displacement
    double areaphase; 		// actual phase modal %
	double EquilPhaseFrac; 	// target phase modal %
	double upper=0.005; 	// upper bound melt variation
	double lower=0.001; 	// lowe bound melt variation
	double fudge; 			// damping factor for out of equilibrium phases
	double temp,press; 		// current T & P
	int MaxPhase=85; 		// ID of highest phase used (assumes all phases betweeen 81 & MaxPhase are present, will have to be redone!)
	double phasephaseEnergy,hostphaseEnergy, hosthostEnergy; // (surface energy pairs three types of boundary)
	int no_modes=MaxPhase-80; // number of phases treated
	double ratio; 			// position of T vs bounding Perplex results
	vector<int> seq; 		// array used for randomising node selection


    if (ElleCount()==0) ElleAddDoubles(); // add double nodes if gaps too large
    if (ElleDisplay()) EllePlotRegions(ElleCount()); // display input file

    incroffset = ElleSwitchdistance() * 0.01; // trial node displacement distance

    max = ElleMaxNodes();

	for (j=0;j<max;j++)
		seq.push_back(j);    // make list of node numbers to be randomised later on

    for (k=0;k<max;k++) { // checks to make sure double node & triple node spacing OK
         if (ElleNodeIsActive(k)) {
            if (ElleNodeIsDouble(k))
               ElleCheckDoubleJ(k);
            else if (ElleNodeIsTriple(k))
               ElleCheckTripleJ(k);
         }
    }

	modal_stats=fopen("modal_stats.txt","w"); // file to store output modal % stats (reset each time process is run)

    for (i=0;i<EllemaxStages();i++) // loop through for all time steps
	{
		temp=ElleTemperature(); // get current T conditions from Elle file
		for(j=0;j<no_pts;j++) // calc modal % for current PT conditions using linear interpolation
		{
			if((modal[j][0] > temp && modal[j+1][0] <= temp )  || (modal[j][0] < temp && modal[j+1][0] >= temp ) )
			{
				ratio=(modal[j][0]-temp)/(modal[j][0]-modal[j+1][0]);

				for(k=1;k<=no_modes;k++)
					modal_out[k]=(ratio*(modal[j+1][k]-modal[j][k]))+modal[j][k];
				break;
			}
		}

		fprintf(modal_stats,"%lf",temp);

		for(PhaseId=81;PhaseId<=MaxPhase;PhaseId++) // for each phase, update gb positions
		{

			EquilPhaseFrac = modal_out[min_names[PhaseId-80]]/100; //convert predicted interpolated modal % to decimal form
			areaphase = Calcareaphase(PhaseId); // get current area of specified phase

			phasephaseEnergy=surface_energy[PhaseId-81][0]; // read surface energies from table, this will need to read an external database
			hostphaseEnergy=surface_energy[PhaseId-81][1];
			hosthostEnergy=surface_energy[PhaseId-81][2];

			if((fabs(areaphase-EquilPhaseFrac))/EquilPhaseFrac>upper) // update damping factor up or down
				CurrData[Magic]=CurrData[Magic]*1.0075;
			else if((fabs(areaphase-EquilPhaseFrac))/EquilPhaseFrac<lower)
				CurrData[Magic]=CurrData[Magic]/1.0025;


			if((fabs(areaphase-EquilPhaseFrac))/EquilPhaseFrac>upper)
				printf("T phase i areaphase equil %%diseq fudge : %lf %d %d %lf %lf %lf %lf\n", temp, PhaseId,i,areaphase,EquilPhaseFrac,
            	100*(fabs(areaphase-EquilPhaseFrac)/EquilPhaseFrac),
            	CurrData[Magic]); // output to screen to check phase area changes

			fprintf(modal_stats," %lf %lf",EquilPhaseFrac,areaphase); // dump predicted vs actual data for this phase

        	max = ElleMaxNodes();


			random_shuffle(seq.begin(),seq.end()); // randomise order of node selection

			for (n=0;n<max;n++)
			{
            	k=seq[n];

                if (ElleNodeIsActive(k)) // ie node exists
				{
                    Phase_MoveNode(k,&incr,&same,areaphase,CurrData[Mobility],CurrData[Magic],
                                   EquilPhaseFrac,phasephaseEnergy,hostphaseEnergy,hosthostEnergy); // calc node movement

					ElleCrossingsCheck(k,&incr); // actually move the node, and check for crossings
                 }
            }

        	for (j=0;j<max;j++) { // checks to make sure double node & triple node spacing OK
            	if (ElleNodeIsActive(j)){
                	if (ElleNodeIsDouble(j)) ElleCheckDoubleJ(j);
                	else if (ElleNodeIsTriple(j)) ElleCheckTripleJ(j);

            	}
        	}

    	}

		fprintf(modal_stats,"\n");
		fflush(modal_stats);

		// update T & P by increments defined from command line

		temp=ElleTemperature();
		temp+=CurrData[IncT];
		ElleSetTemperature(temp);

		press=EllePressure();
		press+=CurrData[IncP];
		ElleSetPressure(press);

		ElleUpdate();
	}
	fclose(modal_stats);
}

int Phase_MoveNode(int n, Coords *movedir, int *same, double areaphase,
                 double factor, double fudge,double areaequil,
                 double energyphasephase, double energyhostphase,
                 double energyhosthost)
{
    register int    i,l;
    double     tryenergy,newenergy,dx,dy,dx2,dy2,dEdX,mobility,startenergy;
    double     temp,velocity,distance;
    double   randangle,angle;
    Coords   xy, newxy,incr,prev,new_incr;
    ERegion  rgn[3];
    double   new_area[3],area[3],old_area[3] ;
    double   old_density[3], new_density[3], density[3], energy[3] ;
    double   surfaceenergy=0,bodyenergy=0;
    double   startsurfaceenergy=0,startbodyenergy=0;
    int      themin,mintype[3];
    double   minmobfactor;


	// so far there are no body energies in this model, so these values are always zero

    ElleNodePosition(n,&xy);
    newxy = xy;

    ElleRegions(n,rgn);

    mobilityofboundary=factor; // user data 1

    randangle=(rand()/(double)RAND_MAX/2.0)*M_PI;
    dx=cos(randangle)*incroffset;
    dy=sin(randangle)*incroffset;
    //printf("rand=%lf\tdx=%le\ty=%le\n",randangle*180.0/M_PI,dx,dy);

	double dblareaphase = Calcareaphase(PhaseId);
    GetSurfaceNodeEnergy(n,PhaseId,&newxy,&startsurfaceenergy,dblareaphase,
                         fudge,areaequil,energyphasephase,energyhostphase,
                         energyhosthost);
    //printf("startbody=%le\tstartsurface=%le\n",startbodyenergy,startsurfaceenergy);
    startenergy=startbodyenergy+startsurfaceenergy;
    movedir->x=0;
    movedir->y=0;
    newenergy = 0;

    ElleNodePosition(n,&xy);
    ElleNodePrevPosition(n,&prev);
    for(i=0;i<4;i++)
    {
        angle=i*M_PI/2;
        rotate_coords(dx, dy, 0.0, 0.0, &dx2, &dy2, angle);

        newxy.x = xy.x+dx2;
        newxy.y = xy.y+dy2;
        ElleSetPosition(n,&newxy);

	    dblareaphase = Calcareaphase(PhaseId);
        GetSurfaceNodeEnergy(n,PhaseId,&newxy,&surfaceenergy,dblareaphase,
                             fudge,areaequil,energyphasephase,energyhostphase,
                             energyhosthost);

        tryenergy=(bodyenergy-startbodyenergy)+(surfaceenergy-startsurfaceenergy);

        if (tryenergy<newenergy)
        {
            newenergy=tryenergy;
            new_incr.x=movedir->x=dx2;
            new_incr.y=movedir->y=dy2;
        }
    }
    ElleSetPosition(n,&xy);
    ElleSetPrevPosition(n,&prev);

    dEdX=-newenergy;

    //printf("dE=%le\n",dEdX);
    if (dEdX>0.0)
    {
        temp=ElleTemperature();
        mobility=mobilityofboundary*ElleSpeedup();
        velocity=mobility*dEdX;
        distance=velocity*truetimestep;

        if (distance>1.0)
        {
            distance=1.0;
            printf("distance=%le\n",distance);
        }
        movedir->x *= distance/(incroffset*lengthscale);
           movedir->y *= distance/(incroffset*lengthscale);
        new_incr.x *=distance/(incroffset*lengthscale);
        new_incr.y *=distance/(incroffset*lengthscale);
        //printf("io=%le temp=%le dE=%le mob=%le velocity=%le distance=%le x=%le y=%le\n",
            //incroffset,temp,dEdX,mobility,velocity,distance,new_incr.x,new_incr.y);
        newxy.x = xy.x+new_incr.x;
        newxy.y = xy.y+new_incr.y;
        ElleSetPosition(n,&newxy);
        if(ElleFlynnAttributeActive(DISLOCDEN))
        {
            for (l=0;l<3;l++) {
                if (rgn[l]!=NO_NB) {
                    new_area[l] = fabs((double)ElleRegionArea (rgn[l]));
                    if (new_area[l]>old_area[l]){
								new_density[l] = (old_area[l]/new_area[l])*old_density[l];
                         //printf("area ratio= %e\n",old_area[l]/new_area[l]);
                        ElleSetFlynnRealAttribute(rgn[l],new_density[l],DISLOCDEN);
                        //^****************************EEEEEEEEEKKKKKKK!!!!!!!!!
                    }
                }
            }
        }
        ElleSetPosition(n,&xy);
        ElleSetPrevPosition(n,&prev);
    }
    else
    {
        movedir->x = 0;
        movedir->y = 0;
    }
}

int GetBodyNodeEnergy(int n, double *total_energy) // not used in this process (yet)
{
    int        i,nb[3];
    ERegion rgn[3];
    double  sum_energy;
    double   area[3],energy[3];
    double density[3];
    UserData udata;
    ElleRegions(n,rgn);

    sum_energy =0;

    for (i=0;i<3;i++) {
        if (rgn[i]!=NO_NB) {
            ElleGetFlynnRealAttribute(rgn[i],&density[i],DISLOCDEN);
            area[i] = fabs(ElleRegionArea (rgn[i]))*(lengthscale*lengthscale);
            energy[i] = energyofdislocations*area[i]*density[i]*dislocationdensityscaling;
            sum_energy += energy[i];
            //printf("area=%le\tdensity=%le\tenergy=%le\n",area[i],density[i],energy[i]);

        }
    }

    *total_energy=sum_energy;
    //*total_energy=0.0; //  <*****************************EEEEEEEEEEEEEEKKKKKK

    return(0);
}

int GetCSLFactor(Coords_3D xyz[3], float *factor)
{
    double misorientation;
    double tmp;
    double dotprod;

    dotprod=(xyz[0].x*xyz[1].x)+(xyz[0].y*xyz[1].y)+(xyz[0].z*xyz[1].z);

    if(fabs(dotprod-1.0) < 0.0001)
    	misorientation=0.0;
    else
    	misorientation = fabs(acos(dotprod)*RTOD);

    if(misorientation>90.0)
    	misorientation=180-misorientation;

    tmp=1.0/((misorientation/90.0)+0.1);

    tmp=(10.0-tmp)/10.0;

    //*factor=(float)tmp;  // gbm version (supresses sub-gbm)
    //*factor=1.0-tmp; // sub gbm version (supresses gbm)


    // CSL function produces 0% at 0 degrees, 95% at 10 degrees
    // and 99% at 20 degrees c axis misorientation

	*factor=1.0; // normal mode, no CSL
}


#if XY
//function shifted to basecode
int ElleGetFlynnEulerCAxis(int flynn_no, Coords_3D *dircos);
int ElleGetFlynnEulerCAxis(int flynn_no, Coords_3D *dircos)
{
    double alpha,beta,gamma;

    ElleGetFlynnEuler3(flynn_no,&alpha,&beta,&gamma);

	alpha=alpha*DTOR;
	beta=beta*DTOR;

/*
 * Euler angles in ZYZ convention
    dircos->x=sin(beta)*cos(alpha);
    dircos->y=-sin(beta)*sin(alpha);
    dircos->z=cos(beta);
*/

    //converts Euler angles into direction cosines of c-axis
/*
 * Euler angles in ZXZ convention
*/
    dircos->x=sin(beta)*sin(alpha);
    dircos->y=-sin(beta)*cos(alpha);
    dircos->z=cos(beta);
}
#endif

int GetSurfaceNodeEnergy(int n,int meltmineral,Coords *xy,double *energy,
                        double areaphase, double fudge, double areaequil,
                        double energyphasephase, double energyhostphase,
                        double energyhosthost)
{
    int        i,nb[3];
    Coords  xynb;
    double    cosalpha,alpha,E;
    double  len;
    int        G0,G1;
    ERegion  rgn[3], rgn2,rgn3;
    Coords_3D xyz[3];
    float csl_factor;
    int val2,val3;
    double energyofsurface;

    for (i=0;i<3;i++)
    	xyz[i].x = xyz[i].y = xyz[i].z = 0;

    ElleRegions(n,rgn);
    ElleNeighbourNodes(n,nb);
    E = 0.0;

    for (i=0;i<3;i++) {
        if (nb[i]!=NO_NB) {
            ElleRelPosition(xy,nb[i],&xynb,&len);
            ElleNeighbourRegion(nb[i],n,&rgn2); // get one neighbouring flynn
            ElleNeighbourRegion(n,nb[i],&rgn3); // get other neighbouring flynn


	    	// meltmineral=1, melt
	    	// meltmineral=0, xl

           // get code of one neighbouring flynn
	    	ElleGetFlynnIntAttribute(rgn2, &val2, MINERAL);
           // get code of other neighbouring flynn
	    	ElleGetFlynnIntAttribute(rgn3, &val3, MINERAL);

    	    energyofsurface=energyhostphase;      // interface energy of xl-melt

	    	if(val2==val3)
	    	{
        		if ( val2==0)
					energyofsurface=energyphasephase; // interface energy of xl-xl
        		else
					energyofsurface=energyhosthost; //interface energy of melt-melt
    	    }

            if(ElleFlynnHasAttribute(rgn2,CAXIS) && ElleFlynnHasAttribute(rgn3,CAXIS) ) // calculate misorientations between grains (currently not used)
            {
        		ElleGetFlynnCAxis(rgn3,&xyz[0]);
        		ElleGetFlynnCAxis(rgn2,&xyz[1]);
        		GetCSLFactor(xyz,&csl_factor);
            }
            else if(ElleFlynnHasAttribute(rgn2,CAXIS))
            {
        		ElleGetFlynnCAxis(rgn2,&xyz[0]);
        		csl_factor=1.0;
            }
            else if(ElleFlynnHasAttribute(rgn3,CAXIS) )
            {
        		ElleGetFlynnCAxis(rgn3,&xyz[0]);
        		csl_factor=1.0;
            }
            else if(ElleFlynnHasAttribute(rgn2,EULER_3)&& ElleFlynnHasAttribute(rgn3,EULER_3))
            {
        		ElleGetFlynnEulerCAxis(rgn3,&xyz[0]);
        		ElleGetFlynnEulerCAxis(rgn2,&xyz[1]);
        		GetCSLFactor(xyz,&csl_factor);
            }
            else if(ElleFlynnHasAttribute(rgn2,EULER_3))
            {
        		ElleGetFlynnEulerCAxis(rgn2,&xyz[0]);
        		csl_factor=1.0;
            }
            else if(ElleFlynnHasAttribute(rgn3,EULER_3) )
            {
        		ElleGetFlynnEulerCAxis(rgn3,&xyz[0]);
        		csl_factor=1.0;
            }
		    else
        		csl_factor=1.0;


            csl_factor=1.0; //  <**************************EEEEEEEEEEEEEEKKKKKK

            if (len==0.0) len = 1.0;

            cosalpha = ((double)xynb.x * xyz[0].x + xynb.y * xyz[0].y)/len;

            /*
             * make sure cosalpha is in range -1,1 for acos
             */
            if (cosalpha > 1.0) cosalpha = 1.0;
            if (cosalpha < -1.0) cosalpha = -1.0;
            if (cosalpha >= 0.0) alpha = acos(cosalpha);
            else                 alpha = acos(-cosalpha);
            E += len*GetAngleEnergy(alpha)*csl_factor*energyofsurface*lengthscale;

	    	//if(val2 != PhaseId || val3 != PhaseId)
			//	printf("alpha xyz[0].x xyz[0].y :%lf %lf %lf\n",alpha*180/3.1415927,xyz[0].x,xyz[0].y );

            cosalpha = ((double)xynb.x * xyz[1].x + xynb.y * xyz[1].y)/len;

        	/*
             * make sure cosalpha is in range -1,1 for acos
             */
            if (cosalpha > 1.0) cosalpha = 1.0;
            if (cosalpha < -1.0) cosalpha = -1.0;
            if (cosalpha >= 0.0) alpha = acos(cosalpha);
            else                 alpha = acos(-cosalpha);

            E += len*GetAngleEnergy(alpha)*csl_factor*energyofsurface*lengthscale;

        }
    }

	E += fabs((areaphase - areaequil)/areaequil) * fudge; // area balancing fudge to ensure global melt fraction is roughly stable


    *energy=(float)E;

    return(0);
}

double GetAngleEnergy(double angle)
{
    int angdeg;
    double    C,energy,sinangle;

    angdeg = (int)(angle*RTOD + 0.5);
    energy = ElleEnergyLUTValue(angdeg);
    return(energy);
}

double Calcareaphase(int meltmineral)
{
    int i;
    int max = ElleMaxFlynns();
    int val;
	double area=0;

    if (ElleFlynnAttributeActive(MINERAL)) {
        for (i=0;i<max;i++)  {
            if (ElleFlynnIsActive(i)) {
                ElleGetFlynnIntAttribute(i,&val,MINERAL);
                if (val==meltmineral) area+=ElleFlynnArea(i);
            }
        }
    }
    return(area);
}

void ElleNodeSameAttrib(int node, int *same, double *sameval, int attrib_id)
{
    int i,j;
    double theval=0;
    double dval[3];
    ERegion rgn[3];


    *same=1;
    *sameval=0;
    if (ElleFlynnAttributeActive(attrib_id))
	{
        ElleRegions(node,rgn);

        for (i=0;i<3;i++)
		{
            dval[i]=0;
            if (rgn[i]!=NO_NB)
			{
                ElleGetFlynnRealAttribute(rgn[i], &dval[i], attrib_id);
                theval=dval[i];
            }
        }

        for(i=0;i<3;i++)
            if(rgn[i]!=NO_NB && dval[i]!=theval) *same=0;

        if (*same==1) *sameval=theval;
    }
}
