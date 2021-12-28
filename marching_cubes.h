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


	class octonion
	{
	public:
		inline octonion(void)
		{
			for (size_t i = 0; i < 8; i++)
				vertex_data[i] = 0;
		}

		float magnitude(void)
		{
			float result = 0;

			for (size_t i = 0; i < 8; i++)
				result += vertex_data[i] * vertex_data[i];

			return sqrtf(result);
		}


		octonion operator+(const octonion& rhs)
		{
			octonion result;

			for (size_t i = 0; i < 8; i++)
				result.vertex_data[i] = vertex_data[i] + rhs.vertex_data[i];

			return result;
		}

		inline octonion(
			float src_r,
			float src_i,
			float src_j,
			float src_k,
			float src_u1,
			float src_i1,
			float src_j1,
			float src_k1)
		{
			vertex_data[0] = src_r;
			vertex_data[1] = src_i;
			vertex_data[2] = src_j;
			vertex_data[3] = src_k;
			vertex_data[4] = src_u1;
			vertex_data[5] = src_i1;
			vertex_data[6] = src_j1;
			vertex_data[7] = src_k1;
		}

		float vertex_data[8];
	};

	octonion trad_mul(const octonion& qA, const octonion& qB);

	vertex_3 vertex_interp(const float isovalue, vertex_3 p1, vertex_3 p2, float valp1, float valp2);
	short unsigned int tesselate_grid_cube(octonion C, float z_w, const float isovalue, const float upper_threshold, const float lower_threshold, short unsigned int max_iterations, const grid_cube& grid, triangle* const triangles);
	void tesselate_adjacent_xy_plane_pair(octonion C, float z_w, const float isovalue, const float upper_threshold, const float lower_threshold, short unsigned int max_iterations, vector<float> xyplane0, vector<float> xyplane1, const size_t z, vector<triangle> &triangles, const float x_grid_min, const float x_grid_max, const size_t x_res, const float y_grid_min, const float y_grid_max, const size_t y_res, const float z_grid_min, const float z_grid_max, const size_t z_res);


	octonion sin(const octonion& in);

	octonion exp(const octonion& in);

	octonion ln(const octonion& in);

	octonion mul(const octonion& in_a, const octonion& in_b);





	quintonion sin(const quintonion& in);

	quintonion exp(const quintonion& in);

	quintonion ln(const quintonion& in);

	quaternion traditional_mul(const quaternion& in_a, const quaternion& in_b);

	quintonion mul(const quintonion& in_a, const quintonion& in_b);

	quintonion conj_number_type(quintonion& in);

	quintonion pow_number_type(quintonion& in, float exponent);




	inline float iterate(
		octonion Z,
		octonion C,
		float z_w,
		const short unsigned int max_iterations,
		const float threshold);

	vertex_3 vertex_interp_refine(
		octonion C,
		float z_w,
		float isovalue,
		float upper_threshold,
		float lower_threshold,
		short unsigned int max_iterations,
		vertex_3 v0, vertex_3 v1,
		float val_v0, float val_v1);

};

#endif
