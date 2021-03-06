// Modified from Paul Bourke, Polygonising a Scalar Field
// Source code by Shawn Halayka	
// Source code is in the public domain


#include "marching_cubes.h"


namespace marching_cubes
{
	int MC_EdgeTable[256]={
	0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
	0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
	0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
	0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
	0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
	0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
	0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
	0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
	0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
	0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
	0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
	0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
	0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
	0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
	0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
	0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0, // mirrors here...
	0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
	0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
	0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
	0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
	0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
	0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
	0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
	0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
	0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
	0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
	0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
	0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
	0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
	0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
	0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
	0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0   };

	int MC_TriTable[256][16] =
	{{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
	{3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
	{3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
	{3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
	{9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
	{9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
	{2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
	{8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
	{9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
	{4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
	{3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
	{1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
	{4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
	{4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
	{9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
	{5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
	{2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
	{9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
	{0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
	{2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
	{10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
	{4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
	{5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
	{5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
	{9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
	{0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
	{1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
	{10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
	{8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
	{2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
	{7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
	{9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
	{2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
	{11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
	{9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
	{5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
	{11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
	{11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
	{1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
	{9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
	{5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
	{2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
	{0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
	{5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
	{6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
	{3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
	{6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
	{5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
	{1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
	{10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
	{6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
	{8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
	{7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
	{3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
	{5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
	{0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
	{9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
	{8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
	{5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
	{0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
	{6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
	{10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
	{10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
	{8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
	{1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
	{3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
	{0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
	{10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
	{3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
	{6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
	{9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
	{8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
	{3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
	{6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
	{0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
	{10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
	{10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
	{2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
	{7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
	{7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
	{2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
	{1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
	{11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
	{8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
	{0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
	{7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
	{10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
	{2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
	{6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
	{7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
	{2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
	{1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
	{10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
	{10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
	{0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
	{7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
	{6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
	{8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
	{9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
	{6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
	{4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
	{10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
	{8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
	{0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
	{1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
	{8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
	{10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
	{4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
	{10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
	{5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
	{11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
	{9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
	{6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
	{7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
	{3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
	{7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
	{9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
	{3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
	{6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
	{9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
	{1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
	{4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
	{7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
	{6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
	{3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
	{0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
	{6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
	{0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
	{11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
	{6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
	{5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
	{9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
	{1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
	{1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
	{10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
	{0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
	{5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
	{10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
	{11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
	{9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
	{7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
	{2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
	{8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
	{9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
	{9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
	{1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
	{9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
	{9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
	{5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
	{0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
	{10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
	{2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
	{0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
	{0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
	{9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
	{5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
	{3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
	{5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
	{8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
	{0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
	{9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
	{1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
	{3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
	{4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
	{9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
	{11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
	{11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
	{2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
	{9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
	{3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
	{1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
	{4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
	{4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
	{0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
	{3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
	{3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
	{0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
	{9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
	{1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}};

};

vertex_3 marching_cubes::vertex_interp(const float isovalue, vertex_3 p1, vertex_3 p2, float valp1, float valp2)
{
	// Sort the vertices so that cracks don't mess up the water-tightness of the mesh.
	// Note: the cracks don't appear if you use doubles instead of floats.
	if (p2 < p1)
	{
		vertex_3 tempv = p2;
		p2 = p1;
		p1 = tempv;

		float tempf = valp2;
		valp2 = valp1;
		valp1 = tempf;
	}

	const float epsilon = 1e-10f;

	if(fabs(isovalue - valp1) < epsilon)
		return(p1);

	if(fabs(isovalue - valp2) < epsilon)
		return(p2);

	if(fabs(valp1 - valp2) < epsilon)
		return(p1);

	float mu = (isovalue - valp1) / (valp2 - valp1);

	return p1 + (p2 - p1)*mu;
}















short unsigned int marching_cubes::tesselate_grid_cube(octonion C, float z_w, const float isovalue, const float upper_threshold, const float lower_threshold, short unsigned int max_iterations, const grid_cube &grid, triangle *const triangles)
{
	short unsigned int cubeindex = 0;

	if(grid.value[0] < isovalue) cubeindex |= 1;
	if(grid.value[1] < isovalue) cubeindex |= 2;
	if(grid.value[2] < isovalue) cubeindex |= 4;
	if(grid.value[3] < isovalue) cubeindex |= 8;
	if(grid.value[4] < isovalue) cubeindex |= 16;
	if(grid.value[5] < isovalue) cubeindex |= 32;
	if(grid.value[6] < isovalue) cubeindex |= 64;
	if(grid.value[7] < isovalue) cubeindex |= 128;

	if(0 == MC_EdgeTable[cubeindex])
		return 0;

	vertex_3 vertlist[12];

	if(MC_EdgeTable[cubeindex] & 1)
		vertlist[0] = vertex_interp_refine(C, z_w, isovalue, upper_threshold, lower_threshold, max_iterations, grid.vertex[0], grid.vertex[1], grid.value[0], grid.value[1]);

	if(MC_EdgeTable[cubeindex] & 2)
		vertlist[1] = vertex_interp_refine(C, z_w, isovalue, upper_threshold, lower_threshold, max_iterations, grid.vertex[1], grid.vertex[2], grid.value[1], grid.value[2]);

	if(MC_EdgeTable[cubeindex] & 4)
		vertlist[2] = vertex_interp_refine(C, z_w, isovalue, upper_threshold, lower_threshold, max_iterations, grid.vertex[2], grid.vertex[3], grid.value[2], grid.value[3]);

	if(MC_EdgeTable[cubeindex] & 8)
		vertlist[3] = vertex_interp_refine(C, z_w, isovalue, upper_threshold, lower_threshold, max_iterations, grid.vertex[3], grid.vertex[0], grid.value[3], grid.value[0]);

	if(MC_EdgeTable[cubeindex] & 16)
		vertlist[4] = vertex_interp_refine(C, z_w, isovalue, upper_threshold, lower_threshold, max_iterations, grid.vertex[4], grid.vertex[5], grid.value[4], grid.value[5]);

	if(MC_EdgeTable[cubeindex] & 32)
		vertlist[5] = vertex_interp_refine(C, z_w, isovalue, upper_threshold, lower_threshold, max_iterations, grid.vertex[5], grid.vertex[6], grid.value[5], grid.value[6]);

	if(MC_EdgeTable[cubeindex] & 64)
		vertlist[6] = vertex_interp_refine(C, z_w, isovalue, upper_threshold, lower_threshold, max_iterations, grid.vertex[6], grid.vertex[7], grid.value[6], grid.value[7]);

	if(MC_EdgeTable[cubeindex] & 128)
		vertlist[7] = vertex_interp_refine(C, z_w, isovalue, upper_threshold, lower_threshold, max_iterations, grid.vertex[7], grid.vertex[4], grid.value[7], grid.value[4]);

	if(MC_EdgeTable[cubeindex] & 256)
		vertlist[8] = vertex_interp_refine(C, z_w, isovalue, upper_threshold, lower_threshold, max_iterations, grid.vertex[0], grid.vertex[4], grid.value[0], grid.value[4]);

	if(MC_EdgeTable[cubeindex] & 512)
		vertlist[9] = vertex_interp_refine(C, z_w, isovalue, upper_threshold, lower_threshold, max_iterations, grid.vertex[1], grid.vertex[5], grid.value[1], grid.value[5]);

	if(MC_EdgeTable[cubeindex] & 1024)
		vertlist[10] = vertex_interp_refine(C, z_w, isovalue, upper_threshold, lower_threshold, max_iterations, grid.vertex[2], grid.vertex[6], grid.value[2], grid.value[6]);

	if(MC_EdgeTable[cubeindex] & 2048)
		vertlist[11] = vertex_interp_refine(C, z_w, isovalue, upper_threshold, lower_threshold, max_iterations, grid.vertex[3], grid.vertex[7], grid.value[3], grid.value[7]);

	short unsigned int ntriang = 0;

	for(short unsigned int i = 0; MC_TriTable[cubeindex][i] != -1; i += 3)
	{
		triangles[ntriang].vertex[0] = vertlist[MC_TriTable[cubeindex][i  ]];
		triangles[ntriang].vertex[1] = vertlist[MC_TriTable[cubeindex][i+1]];
		triangles[ntriang].vertex[2] = vertlist[MC_TriTable[cubeindex][i+2]];
		ntriang++;
	}

	return ntriang;
}

void marching_cubes::tesselate_adjacent_xy_plane_pair(octonion C, float z_w, const float isovalue, const float upper_threshold, const float lower_threshold, short unsigned int max_iterations, vector<float> xyplane0, vector<float> xyplane1, const size_t z, vector<triangle>& triangles, const float x_grid_min, const float x_grid_max, const size_t x_res, const float y_grid_min, const float y_grid_max, const size_t y_res, const float z_grid_min, const float z_grid_max, const size_t z_res)
{
	float avg = (upper_threshold + lower_threshold) / 2.0f;

	for (size_t i = 0; i < xyplane0.size(); i++)
	{
		float a = fabs(avg - xyplane0[i]);
		float b = fabs(avg - upper_threshold);

		float c = a / b;

		xyplane0[i] = c * isovalue;
	}

	for (size_t i = 0; i < xyplane1.size(); i++)
	{
		float a = fabs(avg - xyplane1[i]);
		float b = fabs(avg - upper_threshold);

		float c = a / b;

		xyplane1[i] = c * isovalue;

	}


    const float x_step_size = (x_grid_max - x_grid_min) / (x_res - 1);
    const float y_step_size = (y_grid_max - y_grid_min) / (y_res - 1);
    const float z_step_size = (z_grid_max - z_grid_min) / (z_res - 1);
 
    for(size_t x = 0; x < x_res - 1; x++)
    {
        for(size_t y = 0; y < y_res - 1; y++)
        {
            grid_cube temp_cube;
 
            size_t x_offset = 0;
            size_t y_offset = 0;
            size_t z_offset = 0;
 
            // Setup vertex 0
            x_offset = 0;
            y_offset = 0;
            z_offset = 0;
            temp_cube.vertex[0].x = x_grid_min + ((x+x_offset) * x_step_size);
            temp_cube.vertex[0].y = y_grid_min + ((y+y_offset) * y_step_size);
            temp_cube.vertex[0].z = z_grid_min + ((z+z_offset) * z_step_size);
 
            if(0 == z_offset)
                temp_cube.value[0] = xyplane0[(x + x_offset)*y_res + (y + y_offset)];
            else
                temp_cube.value[0] = xyplane1[(x + x_offset)*y_res + (y + y_offset)];
 
            // Setup vertex 1
            x_offset = 1;
            y_offset = 0;
            z_offset = 0;
            temp_cube.vertex[1].x = x_grid_min + ((x+x_offset) * x_step_size);
            temp_cube.vertex[1].y = y_grid_min + ((y+y_offset) * y_step_size);
            temp_cube.vertex[1].z = z_grid_min + ((z+z_offset) * z_step_size);
 
            if(0 == z_offset)
                temp_cube.value[1] = xyplane0[(x + x_offset)*y_res + (y + y_offset)];
            else
                temp_cube.value[1] = xyplane1[(x + x_offset)*y_res + (y + y_offset)];
 
            // Setup vertex 2
            x_offset = 1;
            y_offset = 0;
            z_offset = 1;
            temp_cube.vertex[2].x = x_grid_min + ((x+x_offset) * x_step_size);
            temp_cube.vertex[2].y = y_grid_min + ((y+y_offset) * y_step_size);
            temp_cube.vertex[2].z = z_grid_min + ((z+z_offset) * z_step_size);
 
            if(0 == z_offset)
                temp_cube.value[2] = xyplane0[(x + x_offset)*y_res + (y + y_offset)];
            else
                temp_cube.value[2] = xyplane1[(x + x_offset)*y_res + (y + y_offset)];
 
            // Setup vertex 3
            x_offset = 0; 
            y_offset = 0;
            z_offset = 1;
            temp_cube.vertex[3].x = x_grid_min + ((x+x_offset) * x_step_size);
            temp_cube.vertex[3].y = y_grid_min + ((y+y_offset) * y_step_size);
            temp_cube.vertex[3].z = z_grid_min + ((z+z_offset) * z_step_size);
 
            if(0 == z_offset)
                temp_cube.value[3] = xyplane0[(x + x_offset)*y_res + (y + y_offset)];
            else
                temp_cube.value[3] = xyplane1[(x + x_offset)*y_res + (y + y_offset)];
 
            // Setup vertex 4
            x_offset = 0;
            y_offset = 1;
            z_offset = 0;
            temp_cube.vertex[4].x = x_grid_min + ((x+x_offset) * x_step_size);
            temp_cube.vertex[4].y = y_grid_min + ((y+y_offset) * y_step_size);
            temp_cube.vertex[4].z = z_grid_min + ((z+z_offset) * z_step_size);
 
            if(0 == z_offset)
                temp_cube.value[4] = xyplane0[(x + x_offset)*y_res + (y + y_offset)];
            else
                temp_cube.value[4] = xyplane1[(x + x_offset)*y_res + (y + y_offset)];
 
            // Setup vertex 5
            x_offset = 1;
            y_offset = 1;
            z_offset = 0;
            temp_cube.vertex[5].x = x_grid_min + ((x+x_offset) * x_step_size);
            temp_cube.vertex[5].y = y_grid_min + ((y+y_offset) * y_step_size);
            temp_cube.vertex[5].z = z_grid_min + ((z+z_offset) * z_step_size);
 
            if(0 == z_offset)
                temp_cube.value[5] = xyplane0[(x + x_offset)*y_res + (y + y_offset)];
            else
                temp_cube.value[5] = xyplane1[(x + x_offset)*y_res + (y + y_offset)];
 
            // Setup vertex 6
            x_offset = 1;
            y_offset = 1;
            z_offset = 1;
            temp_cube.vertex[6].x = x_grid_min + ((x+x_offset) * x_step_size);
            temp_cube.vertex[6].y = y_grid_min + ((y+y_offset) * y_step_size);
            temp_cube.vertex[6].z = z_grid_min + ((z+z_offset) * z_step_size);
 
            if(0 == z_offset)
                temp_cube.value[6] = xyplane0[(x + x_offset)*y_res + (y + y_offset)];
            else
                temp_cube.value[6] = xyplane1[(x + x_offset)*y_res + (y + y_offset)];
 
            // Setup vertex 7
            x_offset = 0;
            y_offset = 1;
            z_offset = 1;
            temp_cube.vertex[7].x = x_grid_min + ((x+x_offset) * x_step_size);
            temp_cube.vertex[7].y = y_grid_min + ((y+y_offset) * y_step_size);
            temp_cube.vertex[7].z = z_grid_min + ((z+z_offset) * z_step_size);
 
            if(0 == z_offset)
                temp_cube.value[7] = xyplane0[(x + x_offset)*y_res + (y + y_offset)];
            else
                temp_cube.value[7] = xyplane1[(x + x_offset)*y_res + (y + y_offset)];
 
            // Generate triangles from cube.
            static triangle temp_triangle_array[5];
 
            short unsigned int number_of_triangles_generated = tesselate_grid_cube(C, z_w, isovalue, upper_threshold, lower_threshold, max_iterations, temp_cube, temp_triangle_array);

            for(short unsigned int i = 0; i < number_of_triangles_generated; i++)
                 triangles.push_back(temp_triangle_array[i]);
        }
    }
}

quintonion marching_cubes::sin(const quintonion& in)
{
	//	float d = in.x * in.x + in.y * in.y + in.z * in.z + in.w * in.w;


	float e = in.vertex_data[1] * in.vertex_data[1] +
		in.vertex_data[2] * in.vertex_data[2] +
		in.vertex_data[3] * in.vertex_data[3] +
		in.vertex_data[4] * in.vertex_data[4];

	//	float l_d = sqrtf(d);
	float l_e = sqrtf(e);

	quintonion out;

	out.vertex_data[0] = std::sin(in.vertex_data[0]) * cosh(l_e);
	out.vertex_data[1] = in.vertex_data[1] / l_e * cos(in.vertex_data[0]) * sinh(l_e);
	out.vertex_data[2] = in.vertex_data[2] / l_e * cos(in.vertex_data[0]) * sinh(l_e);
	out.vertex_data[3] = in.vertex_data[3] / l_e * cos(in.vertex_data[0]) * sinh(l_e);
	out.vertex_data[4] = in.vertex_data[4] / l_e * cos(in.vertex_data[0]) * sinh(l_e);

	return out;
}


marching_cubes::octonion marching_cubes::exp(const octonion& in)
{
	float all_self_dot = 0;
	float imag_self_dot = 0;
	octonion out;

	for (size_t i = 0; i < 8; i++)
		all_self_dot += (in.vertex_data[i] * in.vertex_data[i]);

	for (size_t i = 1; i < 8; i++)
		imag_self_dot += (in.vertex_data[i] * in.vertex_data[i]);

	if (all_self_dot == 0)
	{
		for (size_t i = 0; i < 8; i++)
			out.vertex_data[i] = 0;

		return out;
	}

	const float l_d = sqrtf(all_self_dot);
	const float l_e = sqrtf(imag_self_dot);

	out.vertex_data[0] = std::exp(in.vertex_data[0]) * cos(l_e);

	if (l_e != 0)
	{
		for (size_t i = 1; i < 8; i++)
			out.vertex_data[i] = in.vertex_data[i] / l_e * std::exp(in.vertex_data[0]) * std::sin(l_e);
	}

	return out;
}

marching_cubes::octonion marching_cubes::ln(const octonion& in)
{
	float all_self_dot = 0;
	float imag_self_dot = 0;
	octonion out;

	for (size_t i = 0; i < 8; i++)
		all_self_dot += (in.vertex_data[i] * in.vertex_data[i]);

	for (size_t i = 1; i < 8; i++)
		imag_self_dot += (in.vertex_data[i] * in.vertex_data[i]);

	if (all_self_dot == 0)
	{
		for (size_t i = 0; i < 8; i++)
			out.vertex_data[i] = 0;

		return out;
	}

	const float l_d = sqrtf(all_self_dot);
	const float l_e = sqrtf(imag_self_dot);

	if (in.vertex_data[0] != 0)
	{
		out.vertex_data[0] = log(l_d);
	}

	if (l_e != 0)
	{
		for (size_t i = 1; i < 8; i++)
			out.vertex_data[i] = in.vertex_data[i] / l_e * acos(in.vertex_data[0] / l_d);
	}

	return out;
}


marching_cubes::octonion marching_cubes::mul(const octonion& in_a, const octonion& in_b)
{
	// A*B == exp(ln(A) + ln(B))
	return exp(ln(in_a) + ln(in_b));
}

marching_cubes::octonion marching_cubes::trad_mul(const octonion& qA, const octonion& qB)
{
	octonion out;

	out.vertex_data[0] = qA.vertex_data[0] * qB.vertex_data[0] - qA.vertex_data[1] * qB.vertex_data[1] - qA.vertex_data[2] * qB.vertex_data[2] - qA.vertex_data[3] * qB.vertex_data[3] - qA.vertex_data[4] * qB.vertex_data[4] - qA.vertex_data[5] * qB.vertex_data[5] - qA.vertex_data[6] * qB.vertex_data[6] - qA.vertex_data[7] * qB.vertex_data[7];
	out.vertex_data[1] = qA.vertex_data[0] * qB.vertex_data[1] + qA.vertex_data[1] * qB.vertex_data[0] + qA.vertex_data[2] * qB.vertex_data[3] - qA.vertex_data[3] * qB.vertex_data[2] + qA.vertex_data[4] * qB.vertex_data[5] - qA.vertex_data[5] * qB.vertex_data[4] - qA.vertex_data[6] * qB.vertex_data[7] + qA.vertex_data[7] * qB.vertex_data[6];
	out.vertex_data[2] = qA.vertex_data[0] * qB.vertex_data[2] - qA.vertex_data[1] * qB.vertex_data[3] + qA.vertex_data[2] * qB.vertex_data[0] + qA.vertex_data[3] * qB.vertex_data[1] + qA.vertex_data[4] * qB.vertex_data[6] + qA.vertex_data[5] * qB.vertex_data[7] - qA.vertex_data[6] * qB.vertex_data[4] - qA.vertex_data[7] * qB.vertex_data[5];
	out.vertex_data[3] = qA.vertex_data[0] * qB.vertex_data[3] + qA.vertex_data[1] * qB.vertex_data[2] - qA.vertex_data[2] * qB.vertex_data[1] + qA.vertex_data[3] * qB.vertex_data[0] + qA.vertex_data[4] * qB.vertex_data[7] - qA.vertex_data[5] * qB.vertex_data[6] + qA.vertex_data[6] * qB.vertex_data[5] - qA.vertex_data[7] * qB.vertex_data[4];
	out.vertex_data[4] = qA.vertex_data[0] * qB.vertex_data[4] - qA.vertex_data[1] * qB.vertex_data[5] - qA.vertex_data[2] * qB.vertex_data[6] - qA.vertex_data[3] * qB.vertex_data[7] + qA.vertex_data[4] * qB.vertex_data[0] + qA.vertex_data[5] * qB.vertex_data[1] + qA.vertex_data[6] * qB.vertex_data[2] + qA.vertex_data[7] * qB.vertex_data[3];
	out.vertex_data[5] = qA.vertex_data[0] * qB.vertex_data[5] + qA.vertex_data[1] * qB.vertex_data[4] - qA.vertex_data[2] * qB.vertex_data[7] + qA.vertex_data[3] * qB.vertex_data[6] - qA.vertex_data[4] * qB.vertex_data[1] + qA.vertex_data[5] * qB.vertex_data[0] - qA.vertex_data[6] * qB.vertex_data[3] + qA.vertex_data[7] * qB.vertex_data[2];
	out.vertex_data[6] = qA.vertex_data[0] * qB.vertex_data[6] + qA.vertex_data[1] * qB.vertex_data[7] + qA.vertex_data[2] * qB.vertex_data[4] - qA.vertex_data[3] * qB.vertex_data[5] - qA.vertex_data[4] * qB.vertex_data[2] + qA.vertex_data[5] * qB.vertex_data[3] + qA.vertex_data[6] * qB.vertex_data[0] - qA.vertex_data[7] * qB.vertex_data[1];
	out.vertex_data[7] = qA.vertex_data[0] * qB.vertex_data[7] - qA.vertex_data[1] * qB.vertex_data[6] + qA.vertex_data[2] * qB.vertex_data[5] + qA.vertex_data[3] * qB.vertex_data[4] - qA.vertex_data[4] * qB.vertex_data[3] - qA.vertex_data[5] * qB.vertex_data[2] + qA.vertex_data[6] * qB.vertex_data[1] + qA.vertex_data[7] * qB.vertex_data[0];

	//cout << qA.r * qB.k1 << " " << qA.i * qB.vertex_data[2]1 << " " << qA.vertex_data[2] * qB.i1 << " " << qA.k * qB.u1 << " " << qA.u1 * qB.k << " " << qA.i1 * qB.vertex_data[2] << " " << qA.vertex_data[2]1 * qB.i << " " << qA.k1 * qB.r << endl;


	return out;
}


marching_cubes::octonion marching_cubes::sin(const octonion& in)
{
	//	float d = in.x * in.x + in.y * in.y + in.z * in.z + in.w * in.w;


	float e = in.vertex_data[1] * in.vertex_data[1] +
		in.vertex_data[2] * in.vertex_data[2] +
		in.vertex_data[3] * in.vertex_data[3] +
		in.vertex_data[4] * in.vertex_data[4] +
		in.vertex_data[5] * in.vertex_data[5] +
		in.vertex_data[6] * in.vertex_data[6] +
		in.vertex_data[7] * in.vertex_data[7];


	//	float l_d = sqrtf(d);
	float l_e = sqrtf(e);

	octonion out;

	out.vertex_data[0] = std::sin(in.vertex_data[0]) * cosh(l_e);

	if (l_e != 0)
	{
		out.vertex_data[1] = in.vertex_data[1] / l_e * cos(in.vertex_data[0]) * sinh(l_e);
		out.vertex_data[2] = in.vertex_data[2] / l_e * cos(in.vertex_data[0]) * sinh(l_e);
		out.vertex_data[3] = in.vertex_data[3] / l_e * cos(in.vertex_data[0]) * sinh(l_e);
		out.vertex_data[4] = in.vertex_data[4] / l_e * cos(in.vertex_data[0]) * sinh(l_e);
		out.vertex_data[5] = in.vertex_data[5] / l_e * cos(in.vertex_data[0]) * sinh(l_e);
		out.vertex_data[6] = in.vertex_data[6] / l_e * cos(in.vertex_data[0]) * sinh(l_e);
		out.vertex_data[7] = in.vertex_data[7] / l_e * cos(in.vertex_data[0]) * sinh(l_e);
	}

	return out;
}



quintonion marching_cubes::exp(const quintonion& in)
{
	//	float d = in.x * in.x + in.y * in.y + in.z * in.z + in.w * in.w;
	float d = in.vertex_data[0] * in.vertex_data[0] +
		in.vertex_data[1] * in.vertex_data[1] +
		in.vertex_data[2] * in.vertex_data[2] +
		in.vertex_data[3] * in.vertex_data[3] +
		in.vertex_data[4] * in.vertex_data[4];

	float e = in.vertex_data[1] * in.vertex_data[1] +
		in.vertex_data[2] * in.vertex_data[2] +
		in.vertex_data[3] * in.vertex_data[3] +
		in.vertex_data[4] * in.vertex_data[4];

	float l_d = sqrtf(d);
	float l_e = sqrtf(e);

	quintonion out;

//	if (in.vertex_data[0] != 0)
//	{
		out.vertex_data[0] = std::exp(in.vertex_data[0]) * cos(l_e);
//	}

	if (l_e != 0)
	{
		out.vertex_data[1] = in.vertex_data[1] / l_e * std::exp(in.vertex_data[0]) * std::sin(l_e);
		out.vertex_data[2] = in.vertex_data[2] / l_e * std::exp(in.vertex_data[0]) * std::sin(l_e);
		out.vertex_data[3] = in.vertex_data[3] / l_e * std::exp(in.vertex_data[0]) * std::sin(l_e);
		out.vertex_data[4] = in.vertex_data[4] / l_e * std::exp(in.vertex_data[0]) * std::sin(l_e);
	}

	return out;
}

quintonion marching_cubes::ln(const quintonion& in)
{
	float d = in.vertex_data[0] * in.vertex_data[0] +
		in.vertex_data[1] * in.vertex_data[1] +
		in.vertex_data[2] * in.vertex_data[2] +
		in.vertex_data[3] * in.vertex_data[3] +
		in.vertex_data[4] * in.vertex_data[4];

	float e = in.vertex_data[1] * in.vertex_data[1] +
		in.vertex_data[2] * in.vertex_data[2] +
		in.vertex_data[3] * in.vertex_data[3] +
		in.vertex_data[4] * in.vertex_data[4];

	float l_d = sqrtf(d);
	float l_e = sqrtf(e);

	quintonion out;

	if (in.vertex_data[0] != 0)
	{
		out.vertex_data[0] = log(l_d);
	}

	if (l_e != 0)
	{
		out.vertex_data[1] = in.vertex_data[1] / l_e * acos(in.vertex_data[0] / l_d);
		out.vertex_data[2] = in.vertex_data[2] / l_e * acos(in.vertex_data[0] / l_d);
		out.vertex_data[3] = in.vertex_data[3] / l_e * acos(in.vertex_data[0] / l_d);
		out.vertex_data[4] = in.vertex_data[4] / l_e * acos(in.vertex_data[0] / l_d);
	}

	return out;
}

quaternion marching_cubes::traditional_mul(const quaternion& in_a, const quaternion& in_b)
{
	quaternion out;

	// perform traditional multiply
	out.x = in_a.x * in_b.x - in_a.y * in_b.y - in_a.z * in_b.z - in_a.w * in_b.w;
	out.y = in_a.x * in_b.y + in_a.y * in_b.x + in_a.z * in_b.w - in_a.w * in_b.z;
	out.z = in_a.x * in_b.z - in_a.y * in_b.w + in_a.z * in_b.x + in_a.w * in_b.y;
	out.w = in_a.x * in_b.w + in_a.y * in_b.z - in_a.z * in_b.y + in_a.w * in_b.x;

	return out;
}

quintonion marching_cubes::mul(const quintonion& in_a, const quintonion& in_b)
{
	// A*B == exp(ln(A) + ln(B))
	quintonion ln_a = ln(in_a);
	quintonion ln_b = ln(in_b);

	quintonion out;

	out.vertex_data[0] = ln_a.vertex_data[0] + ln_b.vertex_data[0];
	out.vertex_data[1] = ln_a.vertex_data[1] + ln_b.vertex_data[1];
	out.vertex_data[2] = ln_a.vertex_data[2] + ln_b.vertex_data[2];
	out.vertex_data[3] = ln_a.vertex_data[3] + ln_b.vertex_data[3];
	out.vertex_data[4] = ln_a.vertex_data[4] + ln_b.vertex_data[4];

	return exp(out);
}


quintonion marching_cubes::conj_number_type(quintonion& in)
{
	quintonion out;

	out.vertex_data[0] = in.vertex_data[0];

	for (size_t i = 1; i < in.vertex_length; i++)
		out.vertex_data[i] = -in.vertex_data[i];

	return out;
}

quintonion marching_cubes::pow_number_type(quintonion& in, float exponent)
{
	const float beta = exponent;

	float d = 0;
	float e = 0;
	quintonion out;

	for (size_t i = 0; i < in.vertex_length; i++)
		d += (in.vertex_data[i] * in.vertex_data[i]);

	for (size_t i = 1; i < in.vertex_length; i++)
		e += (in.vertex_data[i] * in.vertex_data[i]);

	if (d == 0)
	{
		for (size_t i = 0; i < out.vertex_length; i++)
			out.vertex_data[i] = 0;

		return out;
	}

	const float l_d = sqrtf(d);
	const float l_e = sqrtf(e);
	const float d_b2 = powf(d, beta / 2.0f);

	const float theta = beta * acos(in.vertex_data[0] / l_d);

	out.vertex_data[0] = d_b2 * cos(theta);

	if (e != 0)
	{
		for (size_t i = 1; i < out.vertex_length; i++)
			out.vertex_data[i] = (in.vertex_data[i] / l_e) * d_b2 * std::sin(theta);
	}

	return out;
}


vertex_3 marching_cubes::vertex_interp_refine(
	octonion C,
	float z_w,
	float isovalue,
	float upper_threshold,
	float lower_threshold,
	short unsigned int max_iterations,
	vertex_3 v0, vertex_3 v1,
	float val_v0, float val_v1)
{
	if (v1 < v0)
	{
		vertex_3 temp(v0);
		float temp_val = val_v0;

		v0 = v1;
		val_v0 = val_v1;

		v1 = temp;
		val_v1 = temp_val;
	}

	// Start half-way between the vertices.
	vertex_3 result = (v0 + v1) * 0.5f;

	const size_t vertex_refinement_steps = 8;
	const float threshold = isovalue;

	// Refine the result, if need be.
	if (0 < vertex_refinement_steps)
	{
		vertex_3 forward, backward;

		// If p1 is outside of the surface and p2 is inside of the surface ...
		if (val_v0 > val_v1)
		{
			forward = v0;
			backward = v1;
		}
		else
		{
			forward = v1;
			backward = v0;
		}

		for (size_t i = 0; i < vertex_refinement_steps; i++)
		{
			octonion Z;
			Z.vertex_data[0] = result.x;
			Z.vertex_data[1] = result.y;
			Z.vertex_data[2] = result.z;

			float x = iterate(Z, C, z_w, max_iterations, threshold);


			float avg = (upper_threshold + lower_threshold) / 2.0f;

				float a = fabs(avg - x);
				float b = fabs(avg - upper_threshold);

				float c = a / b;

				x = c * threshold;

			// If point is in the set, then move forward by 1/2 of a step, else move backward by 1/2 of a step ...
			if (threshold > x)
			{
				backward = result;
				result += (forward - result) * 0.5f;
			}
			else
			{
				forward = result;
				result += (backward - result) * 0.5f;
			}
		}
	}

	return result;
}

marching_cubes::octonion truncate(const marching_cubes::octonion& in)
{
	marching_cubes::octonion out = in;

	out.vertex_data[5] = 0;
	out.vertex_data[6] = 0;
	out.vertex_data[7] = 0;

	return out;
}


float marching_cubes::iterate(
	octonion Z,
	octonion C,
	float z_w,
	const short unsigned int max_iterations,
	const float threshold)
{
	Z.vertex_data[3] = z_w;
	Z.vertex_data[4] = z_w;
	Z = truncate(Z);

	for (short unsigned int i = 0; i < max_iterations; i++)
	{
		octonion Z_base = Z;

		Z = trad_mul(Z, Z_base);
		Z = truncate(Z);

		Z = mul(Z, Z_base);
		Z = truncate(Z);

		Z = mul(Z, Z_base);
		Z = truncate(Z);

		Z = Z + C;
		Z = truncate(Z);

		if (Z.magnitude() >= threshold)
			break;
	}

	return Z.magnitude();
}