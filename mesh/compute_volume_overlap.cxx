/* Author: Christine Xu, University of North Carolina at Chapel Hill*/

#ifndef COMPUTE_VOLUME_OVERLAP_C
#define COMPUTE_VOLUME_OVERLAP_C

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <memory.h>

#include <compute_volume_overlap.h>

/* --------------------------------------------------------------------------*
 *                    Local data type                                        *
 * --------------------------------------------------------------------------*/

#define MIN(A,B) (((A) < (B)) ? (A) : (B))
#define MAX(A,B) (((A) > (B)) ? (A) : (B))

typedef struct hit {
	double d;
	int sign;
	char boundary;	// 0=none 1=A 2=B
} HIT;

// STATIC hits
hit hits[400][400][100];	// [pix][pix][# of hits]
//2D grid.  each element points to a list of 'hits', distances at which a 
//triangle was found.  Think of french fries.

/* --------------------------------------------------------------------------*
 *                    Local utility functions                                *
 * --------------------------------------------------------------------------*/

//sort function
int compareHitByDistance(const void *arg1, const void *arg2)
{
	HIT *h1 = (HIT*)arg1;
	HIT *h2 = (HIT*)arg2;
	if (h1->d < h2->d)
		return -1;
	if (h1->d > h2->d)
		return  1;
	return 0;
}


//Simple swap function.
void swap(double v[], int ind1, int ind2)
{
	double temp = v[ind1];
	v[ind1] = v[ind2];
	v[ind2] = temp;
}

// Sorts an array of three numbers in ascending order in place.
void sort3( double* a )
{
	double tmp;
	for(int i = 0; i < 2; ++i ) {
		for(int j = i+1; j < 3; ++j ) {
			if( a[i] > a[j] ) {
				tmp = a[j];
				a[j] = a[i];
				a[i] = tmp;
			}
		}
	}
}


/* --------------------------------------------------------------------------*
 *                    Local functions                                        *
 * --------------------------------------------------------------------------*/

// find, in ray (i,j), distances of {A, B, union(A,B), intersection(A,B)};
// RETURN: hit count
//Returns volumes along this grid element, returned in d.
int findDistances(int i, int j, double d[4] )
{
	// no hits on this ray? -> return "no hits"
	if (! hits[i][j][0].boundary)
		return 0;

	HIT *h;	// local hit ptr
	int nHitsA = 0, nHitsB = 0;

	// if A or B have an odd hit count, fail
	h = hits[i][j];
	while (h->boundary)
	{
		if (h->boundary == 1)
			nHitsA++;
		else
			nHitsB++;
		h++;
	}
	if (nHitsA&1  ||  nHitsB&1)
	{
		printf ("WARNING: odd hits on A (%d) or B (%d)\n", nHitsA, nHitsB);
		return 0;	// no hits were processed [LATER: process the non-odd volume]
	}

	// sort hits by distance
	h = hits[i][j];
	int nHits = nHitsA + nHitsB;
	qsort( (void *)h, (size_t)nHits, sizeof(HIT), compareHitByDistance);

	// calc d[]
	// - add this segment's distance, if positive, to proper counters
	// - update to look at next intersection
	int inA = 0, inB = 0;	// is the next segment in vol A? in vol B ?
	double lastDist = -1;	// beginning of last segment
	int addWeight[4][4] =	// [inA*2+inB][a,b,union,intersection]
	{
		{0,0,0,0},	// A:out B:out
		{0,1,1,0},	// A:out B:in
		{1,0,1,0},	// A:in  B:out
		{1,1,1,1}	// A:in  B:in
	};
	d[0] = d[1] = d[2] = d[3] = 0.0;

	h = hits[i][j];
	while (h->boundary)
	{
		// add the distance from the last intersection into the appropriate d[] counters
		if (lastDist != -1)	// there's no segment (to add) before the first one
		{
			for (int c=0; c < 4; c++)
				d[c] += addWeight[inA*2+inB][c] * (h->d - lastDist);
		}
		lastDist = h->d;

		// recalc which surface has just been entered or exited
		if (h->boundary == 1)
			inA = (h->sign < 0);
		else
			inB = (h->sign < 0);
		h++;
	}
	
	return nHits;
}


// add a hit (d,s,b) to first empty slot in ray (i,j)
void addHit(int i, int j, double distance, int sign, int boundary)
{
	// boundary==1 -> use Ahits
	// boundary==2 -> use Bhits
	HIT *h = hits[i][j];

	// find first unused entry
	int nHits = 0;
	while (h->boundary)
	{
		if (nHits++ >= 100) {printf("ERROR: exceeded 100 hits per ray - exiting\n"); exit(1);}
		h++;
	}

	h->boundary = boundary;
	h->d        = distance;
	h->sign     = sign;
}

double Intersect_tri(double e[], double v[], double tuv[], bool Tri1, int trind, double * L1, double * L2, int * T1, int * T2)
//Returns the extent along v from e to hit the triangle (negative if 
//non-intersect). In the end, tuv[0] will contain a scaler for v, and 
//tuv[1] will be a sign, NEGATIVE if the triangle's normal points to e.
//Optimized by assuming orthographic projection.  I return false before 
//much computation if e[0] e[1] is not within the bounding x y of the
//triangle.  This should prevent several dozen multiplications and mallocs.
{
        double edge1[3], edge2[3];
	int i;

	if (Tri1) {
		for (i = 0; i < 3; i ++) {
			edge1[i] = L1[(T1[trind+1]-1)*3+i] - L1[(T1[trind]-1)*3+i]; //t->p1[i] - t->p0[i];
			edge2[i] = L1[(T1[trind+2]-1)*3+i] - L1[(T1[trind]-1)*3+i]; //t->p2[i] - t->p0[i];
		}
	}else {
		for (i = 0; i < 3; i ++) {
			edge1[i] = L2[(T2[trind+1]-1)*3+i] - L2[(T2[trind]-1)*3+i]; //t->p1[i] - t->p0[i];
			edge2[i] = L2[(T2[trind+2]-1)*3+i] - L2[(T2[trind]-1)*3+i]; //t->p2[i] - t->p0[i];
		}
	}


	double tvec[3], pvec[3], qvec[3], det,inv_det;


	//Cross v and edge2.
	pvec[0] = v[1]*edge2[2] - v[2]*edge2[1];
	pvec[1] = v[2]*edge2[0] - v[0]*edge2[2];
	pvec[2] = v[0]*edge2[1] - v[1]*edge2[0];

	//dot edge1 with the cross of v and edge2.  If the dot is close to 
	//zero, the v is in the plane of the triangle.
	det = edge1[0]*pvec[0] + edge1[1]*pvec[1] + edge1[2]*pvec[2];

	if (fabs(det) < 1.0e-10) return -10000.0;

	inv_det = 1.0/det;

	//Distance from vert0 to ray origin.
	if (Tri1) {
		tvec[0] = e[0]- L1[(T1[trind]-1)*3]; //t->p0[0];
		tvec[1] = e[1]- L1[(T1[trind]-1)*3+1]; //t->p0[1];
		tvec[2] = e[2]- L1[(T1[trind]-1)*3+2]; //t->p0[2];
	}else {
		tvec[0] = e[0]- L2[(T2[trind]-1)*3]; //t->p0[0];
		tvec[1] = e[1]- L2[(T2[trind]-1)*3+1]; //t->p0[1];
		tvec[2] = e[2]- L2[(T2[trind]-1)*3+2]; //t->p0[2];
	}

	//Calculate u and test bounds
	tuv[1] = tvec[0]*pvec[0] + tvec[1]*pvec[1] + tvec[2]*pvec[2];
	tuv[1] = tuv[1]*inv_det;
	if ((tuv[1] < -1.0e-10)||(tuv[1] > 1.0)) return -10000.0;

	//Prepare for v parameter.
	qvec[0] = tvec[1]*edge1[2] - tvec[2]*edge1[1];
	qvec[1] = tvec[2]*edge1[0] - tvec[0]*edge1[2];
	qvec[2] = tvec[0]*edge1[1] - tvec[1]*edge1[0];

	//Calculate v and test bounds
	tuv[2] = v[0]*qvec[0] + v[1]*qvec[1] + v[2]*qvec[2];
	tuv[2] = tuv[2]*inv_det;

	if ((tuv[2] < -1.0e-10)||(tuv[1] + tuv[2] > 1.0)) return -10000.0;

	//Calculate t.
	tuv[0] = edge2[0]*qvec[0] + edge2[1]*qvec[1] + edge2[2]*qvec[2];
	tuv[0] = tuv[0]*inv_det;
	
	//We'll now put into tuv[1] the sign of the hit, v dotted with the 
	//normal to the triangle. Reusing variables.

	pvec[0] = edge1[1]*edge2[2] - edge1[2]*edge2[1];
	pvec[1] = edge1[2]*edge2[0] - edge1[0]*edge2[2];
	pvec[2] = edge1[0]*edge2[1] - edge1[1]*edge2[0];

	//pvec is the normal to the triangle. tuv[1] used to be extent along the 
	//u barrycentric vector.
	tuv[1] = v[0]*pvec[0] + v[1]*pvec[1] + v[2]*pvec[2];



	//We also scale the distance t by the length of vector v, so that
	//this t is comparable to those returned by other functions.

	//tuv[0] = tuv[0]/(sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]));

	return tuv[0];

}


//Pengdong Xiao, April 27, 2012: swap T1 and T2
void GetVolumeOverlap(double vols[],double * L1, double * L2, int numVertices1, int numVertices2, int * T1, int * T2, int numTriangles1, int numTriangles2) 
//L1, L2 are vertex lists (x coord, y, z, and then the next point, etc).
//T1 and T2 are index lists for the tileset triangles.  So vert1_index, vert2_index, vert3_index for every triangle 

//In order to approximate the volume overlap of the two polyhedra, 
//we determine a plane underneath the two objects, the size of the
//projected bounding box.  Then we orthographically ray trace the
//objects (bad ray tracing--for every pixel for every triangle.  
//For each ray, we get a line with in and out markers for both 
//volumes and we can take a weighted sum of the distances times
//a differential area (pixel size) as the intersection volume over
//that pixel.  Summed over all pixels, we have the total
//intersection volume. Concurrently compute the volume of the 
//reference to get volume overlap.
{

	int i, j, t, temp, pix = 400;	// CAREFUL! change size of Ahits if you change this!
	double minx = DBL_MAX, miny = DBL_MAX, minz = DBL_MAX;
	double maxx = -DBL_MAX, maxy = -DBL_MAX;  //for the bounding plane.
	double xextent, yextent, xinc, yinc, d[4], dA;
	double sumx = 0, sumy = 0;
	double A_vol = 0, B_vol = 0.0, intersection_vol = 0.0, union_vol = 0;

	double e[3], n[3], ret = 0.0;  //Point and normal of the plane.

	//printf("Computing Volume Overlap Approximation.\n");

	// reset intersection hits
	memset((void*)hits, 0, sizeof(hits));

	temp = 3*numVertices1;
	for (i = 0; i < temp; i = i + 3) {
		if (L1[i] < minx) minx = L1[i];
		if (L1[i] > maxx) maxx = L1[i];
		sumx += L1[i];
		if (L1[i+1] < miny) miny = L1[i+1];
		if (L1[i+1] > maxy) maxy = L1[i+1];
		sumy += L1[i+1];
		if (L1[i+2] < minz) minz = L1[i+2];
	}
	
	temp = 3*numVertices2;
	for (i = 0; i < temp; i = i + 3) {
		if (L2[i] < minx) minx = L2[i];
		if (L2[i] > maxx) maxx = L2[i];
		sumx += L2[i];
		if (L2[i+1] < miny) miny = L2[i+1];
		if (L2[i+1] > maxy) maxy = L2[i+1];
		sumy += L2[i+1];
		if (L2[i+2] < minz) minz = L2[i+2];
	}

	//printf("\tBounding box (x y z): %f %f, %f %f, %f\n", minx,maxx,miny,maxy,minz);

	e[0] = minx; e[1] = miny; e[2] = minz - 5.0;
	n[0] =  0.0; n[1] = 0.0; n[2] = 1.0;
	xextent = maxx-minx; yextent = maxy-miny;
	xinc = xextent/pix; yinc = yextent/pix;
	dA = xinc*yinc;

	//printf("pix = %d, thus column base area = %lf ( %f * %f )\n", pix, dA, xinc, yinc);

	int vert;
	// for each tile, intersect rays that cover it's bounding box with the tile
	for (t = 0; t < numTriangles1; t++)
	{
		// find tile's bb
		vert = T1[3*t + 0]-1;
		double tileMinx =        L1[vert*3 + 0];
		double tileMaxx =        L1[vert*3 + 0];
		double tileMiny =        L1[vert*3 + 1];
		double tileMaxy =        L1[vert*3 + 1];
		
		vert = T1[3*t + 1]-1;
		tileMinx = MIN(tileMinx, L1[vert*3 + 0]);
		tileMaxx = MAX(tileMaxx, L1[vert*3 + 0]);		
		tileMiny = MIN(tileMiny, L1[vert*3 + 1]);
		tileMaxy = MAX(tileMaxy, L1[vert*3 + 1]);
		
		vert = T1[3*t + 2]-1;
		tileMinx = MIN(tileMinx, L1[vert*3 + 0]);
		tileMaxx = MAX(tileMaxx, L1[vert*3 + 0]);		
		tileMiny = MIN(tileMiny, L1[vert*3 + 1]);
		tileMaxy = MAX(tileMaxy, L1[vert*3 + 1]);
		
		// now intersect all rays in bb with the tile and record in Ahits
		double tuv[3], dst;
		int startx = int((tileMinx - minx) / xinc);
		int   endx = int((tileMaxx - minx) / xinc);
		int starty = int((tileMiny - miny) / yinc);
		int   endy = int((tileMaxy - miny) / yinc);
		for (i = startx; i <= endx; i++)
		{
			for (j = starty; j <= endy; j++)
			{
				e[0] = minx + xinc/2.0 + i*xinc;
				e[1] = miny + yinc/2.0 + j*yinc;
				e[2] = minz - 5.0;
				dst = Intersect_tri(e,n,tuv,true, 3*t, L1, L2, T1, T2);
				if (dst > -1.0e3) {	//A real hit, so append to the hitlist.
					int sgn = (int) (tuv[1]/fabs(tuv[1]));
					addHit(i, j, fabs(dst), sgn, 1);
					//numhits1++;
				}
			}							
		}
	}
		
	for (t = 0; t < numTriangles2; t++)
	{
		// find tile's bb
		vert = T2[3*t + 0]-1;
		double tileMinx =        L2[vert*3 + 0];
		double tileMaxx =        L2[vert*3 + 0];
		double tileMiny =        L2[vert*3 + 1];
		double tileMaxy =        L2[vert*3 + 1];
		
		vert = T2[3*t + 1]-1;
		tileMinx = MIN(tileMinx, L2[vert*3 + 0]);
		tileMaxx = MAX(tileMaxx, L2[vert*3 + 0]);		
		tileMiny = MIN(tileMiny, L2[vert*3 + 1]);
		tileMaxy = MAX(tileMaxy, L2[vert*3 + 1]);
		
		vert = T2[3*t + 2]-1;
		tileMinx = MIN(tileMinx, L2[vert*3 + 0]);
		tileMaxx = MAX(tileMaxx, L2[vert*3 + 0]);		
		tileMiny = MIN(tileMiny, L2[vert*3 + 1]);
		tileMaxy = MAX(tileMaxy, L2[vert*3 + 1]);
		
		// now intersect all rays in bb with the tile and record in Ahits
		double tuv[3], dst;
		int startx = int((tileMinx - minx) / xinc);
		int endx = int((tileMaxx - minx) / xinc);
		int starty = int((tileMiny - miny) / yinc);
		int endy = int((tileMaxy - miny) / yinc);
		for (i = startx; i <= endx; i++)
		{
			for (j = starty; j <= endy; j++)
			{
				e[0] = minx + xinc/2.0 + i*xinc;
				e[1] = miny + yinc/2.0 + j*yinc;
				e[2] = minz - 5.0;
				dst = Intersect_tri(e,n,tuv,false, 3*t, L1, L2, T1, T2);
				if (dst > -1.0e3) {	//A real hit, so append to the hitlist.
					int sgn = (int) (tuv[1]/fabs(tuv[1]));
					addHit(i, j, fabs(dst), sgn, 2);
					//numhits1++;
				}
			}							
		}
	}

	for (i = 0; i < pix; i ++) {
		//if (i%10 == 0) {printf(" %d ", i); spittime("");}
		
		for (j = 0; j < pix; j ++) {
			if (findDistances(i,j,d))
			{
				A_vol            += d[0];
				B_vol            += d[1];
				union_vol        += d[2];
				intersection_vol += d[3];
			}
		}
	}

	


	//printf("\n\n");
	//printf(" Volume OVERLAP-A:     %f = %f * base\n", A_vol            * dA, A_vol           );
	//printf(" Volume OVERLAP-B:     %f = %f * base\n", B_vol            * dA, B_vol           );
	//printf(" Volume OVERLAP-UNION: %f = %f * base\n", union_vol        * dA, union_vol       );
	//printf(" Volume OVERLAP-INT:   %f = %f * base\n", intersection_vol * dA, intersection_vol);

	//printf(" INT-UNION:%lf:\n", intersection_vol/union_vol);
	//printf(" INT-AVE:%lf:\n", intersection_vol/((A_vol + B_vol)/2.0));

	//ret = int_vol/B_vol;
	//ret = intersection_vol / B_vol;
	//ret = intersection_vol / ((A_vol + B_vol)/2.0);
			
	
	//return ret;

	vols[0] = A_vol; 
	vols[1] = B_vol;
	vols[2] = union_vol;
	vols[3] = intersection_vol;
}


//coord1 and coord2 are 0, 1, or 2.  This function swaps values in both vertex 
//arrays so that volume overlap projection can be with respect to orthogonal planes.
//In order to preserve the right-hand rule, you should negate one of the dimensions, 
//and I don't think it makes a difference, so I negate the dimension not swapped.
void SwapCoordinates(int coord, double * L1, double * L2, int numVertices1, int numVertices2)
{
	int i; 

	int coord1 = (coord + 1) % 3;
	int coord2 = (coord + 2) % 3;

	double temp = 3*numVertices1;
	for (i = 0; i < temp; i = i + 3)
	{
		swap(L1, (i + coord1), (i + coord2));
		L1[i + coord] = -1.0*L1[i + coord];
	}
	
	temp = 3*numVertices2;
	for (i = 0; i < temp; i = i + 3) 
	{
		swap(L2, (i + coord1), (i + coord2));
		L2[i + coord] = -1.0*L2[i + coord];
	}

}

/* --------------------------------------------------------------------------*
 *                    External functions                                     *
 * --------------------------------------------------------------------------*/

//This is an overarching function that will call GetVolumeOverlap 3 times, one for
//each projection (x,y,z).  That way it's a bit more robust to degenerate or axis-
//aligned byus.
void ComputeRobustVolumeOverlap(double * L1, double * L2, int numVertices1, int numVertices2, int * T1, int * T2, int numTriangles1, int numTriangles2, double dice[], double int_union_ratio[] )
{
	double vols[4];
	double int_avgs[3], int_unions[3];
	double int_avg = 0.0, int_union = 0.0; 

	GetVolumeOverlap(vols, L1, L2, numVertices1, numVertices2, T1, T2, numTriangles1, numTriangles2);
	int_avgs[0]   = vols[3]/((vols[0] + vols[1])/2.0);
	int_unions[0] = vols[3]/vols[2]; 

	SwapCoordinates(1, L1, L2, numVertices1, numVertices2); //becomes (z,y,x)

	GetVolumeOverlap(vols, L1, L2, numVertices1, numVertices2, T1, T2, numTriangles1, numTriangles2);
	int_avgs[1]   = vols[3]/((vols[0] + vols[1])/2.0);
	int_unions[1] = vols[3]/vols[2];

	SwapCoordinates(0, L1, L2, numVertices1, numVertices2); //becomes (z,x,y)

	GetVolumeOverlap(vols, L1, L2, numVertices1, numVertices2, T1, T2, numTriangles1, numTriangles2);
	int_avgs[2]   = vols[3]/((vols[0] + vols[1])/2.0);
	int_unions[2] = vols[3]/vols[2];


	sort3(int_avgs);
	sort3(int_unions);

	int_avg   = (int_avgs[1] + int_avgs[2]) / 2.0;
	int_union = (int_unions[1] + int_unions[2]) / 2.0;
	
	printf("\n\nOVERALL_DICE:%lf\n", int_avg);
	printf("OVERALL_INT/UNION:%lf\n", int_union);

	dice[0] = int_avg;
	int_union_ratio[0] = int_union;
}


#endif
