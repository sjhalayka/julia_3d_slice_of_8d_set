#include "main.h"





int main(void)
{
	vector<pair<float, float>> thresholds;

	vector<vector<triangle>> triangles;

	pair<float, float> f;

	f.first = 4.0f;
	f.second = 0.0f;
	//thresholds.push_back(f);

	//f.first = 0.95f;
	//f.second = 0.55f;
	//thresholds.push_back(f);

	//f.first = 0.5f;
	//f.second = 0.0f;
	//thresholds.push_back(f);

	triangles.resize(thresholds.size());

	const float grid_max = 1.5;
	const float grid_min = -grid_max;
	const size_t res = 100;

	const bool make_border = true;

	const unsigned short int max_iterations = 8;
	const float threshold = 4;

	// When adding a border, use a value that is "much" greater than the threshold.
	const float border_value = 1.0f + threshold;


	vector<float> xyplane0(res * res, 0);
	vector<float> xyplane1(res * res, 0);

	const float step_size = (grid_max - grid_min) / (res - 1);

	octonion C;
	C.vertex_data[0] = 0.3f;
	C.vertex_data[1] = 0.5f;
	C.vertex_data[2] = 0.4f;
	C.vertex_data[3] = 0.2f;
	C.vertex_data[4] = 0.1f;
	C.vertex_data[5] = 0.0f;
	C.vertex_data[6] = 0.0f;
	C.vertex_data[7] = 0.0f;

	octonion Z;

	for (size_t i = 0; i < 3; i++)
		Z.vertex_data[i] = grid_min;

	// Do slice of 8D set
	const float slice_val = 0.45f;

	size_t z = 0;

	// Calculate first xy plane
	for (size_t x = 0; x < res; x++, Z.vertex_data[0] += step_size)
	{
		Z.vertex_data[1] = grid_min;

		for (size_t y = 0; y < res; y++, Z.vertex_data[1] += step_size)
		{
			if (0)//z > res / 2)
				xyplane0[x * res + y] = border_value;
			else
			{
				if (true == make_border && (x == 0 || y == 0 || z == 0 || x == res - 1 || y == res - 1 || z == res - 1))
					xyplane0[x * res + y] = border_value;
				else
					xyplane0[x * res + y] = iterate(Z, C, slice_val, max_iterations, threshold);
			}
		}
	}

	// Prepare for 2nd xy plane
	z++;
	Z.vertex_data[2] += step_size;

	// Calculate 2nd and subsequent xy planes
	for (; z < res; z++, Z.vertex_data[2] += step_size)
	{
		Z.vertex_data[0] = grid_min;

		cout << "Calculating triangles from xy-plane pair " << z << " of " << res - 1 << endl;

		for (size_t x = 0; x < res; x++, Z.vertex_data[0] += step_size)
		{
			Z.vertex_data[1] = grid_min;

			for (size_t y = 0; y < res; y++, Z.vertex_data[1] += step_size)
			{

				if (0)//z > res / 2)
					xyplane1[x * res + y] = border_value;
				else
				{
					if (true == make_border && (x == 0 || y == 0 || z == 0 || x == res - 1 || y == res - 1 || z == res - 1))
						xyplane1[x * res + y] = border_value;
					else
						xyplane1[x * res + y] = iterate(Z, C, slice_val, max_iterations, threshold);
				}
			}
		}



		for (size_t t = 0; t < thresholds.size(); t++)
		{
			// Calculate triangles for the xy-planes corresponding to z - 1 and z by marching cubes
			tesselate_adjacent_xy_plane_pair(
				C, slice_val,
				threshold,
				thresholds[t].first,
				thresholds[t].second,
				max_iterations,
				xyplane0, xyplane1,
				z - 1,
				triangles[t],
				grid_min, grid_max, res,
				grid_min, grid_max, res,
				grid_min, grid_max, res);
		}

		// Swap memory pointers (fast) instead of performing a memory copy (slow)
		xyplane1.swap(xyplane0);
	}

	for (size_t t = 0; t < triangles.size(); t++)
	{
		ostringstream oss;
		oss << t;

		string filename = "out" + oss.str() + ".stl";

		if (0 < triangles[t].size())
			write_triangles_to_binary_stereo_lithography_file(triangles[t], filename.c_str());
	}

	return 0;
}