// Modified from Paul Bourke, Polygonising a Scalar Field
// Source code by Shawn Halayka
// Source code is in the public domain


#ifndef MARCHING_CUBES_H
#define MARCHING_CUBES_H


#include "primitives.h"


#include <iostream>
using std::cout;
using std::endl;

#include <fstream>
using std::ofstream;

#include <iomanip>
using std::ios_base;

#include <vector>
using std::vector;

#include <set>
using std::set;

#include <cmath>
using std::sin;
using std::exp;



namespace marching_cubes
{
	class grid_cube
	{
	public:
		vertex_3 vertex[8];
		float value[8];
	};

//	vertex_3 vertex_interp(const float isovalue, vertex_3 p1, vertex_3 p2, float valp1, float valp2);
	short unsigned int tesselate_grid_cube(quintonion C, float z_w, const float isovalue, const float upper_threshold, const float lower_threshold, short unsigned int max_iterations, const grid_cube& grid, triangle* const triangles);
	void tesselate_adjacent_xy_plane_pair(quintonion C, float z_w, const float isovalue, const float upper_threshold, const float lower_threshold, short unsigned int max_iterations, vector<float> xyplane0, vector<float> xyplane1, const size_t z, vector<triangle> &triangles, const float x_grid_min, const float x_grid_max, const size_t x_res, const float y_grid_min, const float y_grid_max, const size_t y_res, const float z_grid_min, const float z_grid_max, const size_t z_res);


	quintonion sin(const quintonion& in);

	quintonion exp(const quintonion& in);

	quintonion ln(const quintonion& in);

	quaternion traditional_mul(const quaternion& in_a, const quaternion& in_b);

	quintonion mul(const quintonion& in_a, const quintonion& in_b);


	quintonion conj_number_type(quintonion& in);

	quintonion pow_number_type(quintonion& in, float exponent);




	inline float iterate(
		quintonion Z,
		quintonion C,
		float z_w,
		const short unsigned int max_iterations,
		const float threshold)
	{
		Z.vertex_data[3] = z_w;
		Z.vertex_data[4] = z_w;

		for (short unsigned int i = 0; i < max_iterations; i++)
		{
			//quintonion Z_orig = Z;

			quintonion Z_base = Z;
			Z = mul(Z, Z_base);
			Z = mul(Z, Z_base);
			Z = mul(Z, Z_base);

			Z = Z + C;


		//	Z = pow_number_type(Z_orig, 2.0) + C;




		//	Z = sin(Z) + mul(sin(Z), C);


	/*		quaternion qc;
			qc.x = C.vertex_data[0];
			qc.y = C.vertex_data[1];
			qc.z = C.vertex_data[2];
			qc.w = C.vertex_data[3];

			quintonion s = sin(Z);

			quaternion qs;
			qs.x = s.vertex_data[0];
			qs.y = s.vertex_data[1];
			qs.z = s.vertex_data[2];
			qs.w = s.vertex_data[3];

			quaternion f = traditional_mul(qc, qs);

			f.x += qs.x;
			f.y += qs.y;
			f.z += qs.z;
			f.w += qs.w;

			Z.vertex_data[0] = f.x;
			Z.vertex_data[1] = f.y;
			Z.vertex_data[2] = f.z;
			Z.vertex_data[3] = f.w;
			Z.vertex_data[4] = 0;*/


			if (Z.magnitude() >= threshold)
				break;
		}

		return Z.magnitude();
	}

	vertex_3 vertex_interp_refine(
		quintonion C,
		float z_w,
		float isovalue,
		float upper_threshold,
		float lower_threshold,
		short unsigned int max_iterations,
		vertex_3 v0, vertex_3 v1,
		float val_v0, float val_v1);

};

#endif
