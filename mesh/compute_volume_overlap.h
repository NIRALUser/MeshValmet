/* Author: Christine Xu, University of North Carolina at Chapel Hill*/

#ifndef _COMPUTE_VOLUME_OVERLAP_PROTO
#define _COMPUTE_VOLUME_OVERLAP_PROTO



#ifdef __cplusplus
#define BEGIN_DECL extern "C" {
#define END_DECL }
#else
#define BEGIN_DECL
#define END_DECL
#endif

BEGIN_DECL
#undef BEGIN_DECL

//This is an overarching function that will call GetVolumeOverlap 3 times, one for
//each projection (x,y,z).  That way it's a bit more robust to degenerate or axis-
//aligned byus.
void ComputeRobustVolumeOverlap(double * L1, double * L2, int numVertices1, int numVertices2, int * T1, int * T2, int numTriangles1, int numTriangles2, double dice[], double int_union_ratio[] );

END_DECL
#undef END_DECL

#endif /*_COMPUTE_VOLUME_OVERLAP_PROTO*/
