#ifndef GEOM3D_H
#define GEOM3D_H

#define GEOM3D_GENERATE_INDICES 1
#define GEOM3D_GENERATE_POSITIONS 2
#define GEOM3D_GENERATE_NORMALS 4
#define GEOM3D_GENERATE_UVS 8
#define GEOM3D_GENERATE_ALL (GEOM3D_GENERATE_INDICES | GEOM3D_GENERATE_POSITIONS | GEOM3D_GENERATE_NORMALS | GEOM3D_GENERATE_UVS)
#define GEOM3D_GENERATE_NO_UVS (GEOM3D_GENERATE_INDICES | GEOM3D_GENERATE_POSITIONS | GEOM3D_GENERATE_NORMALS)
#define GEOM3D_GENERATE_BASIC (GEOM3D_GENERATE_INDICES | GEOM3D_GENERATE_POSITIONS)

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

	extern int box_buffer_geometry(
		uint8_t options,
		uint32_t **indices, uint32_t *indices_count,
		float **positions, uint32_t *positions_count,
		float **normals, uint32_t *normals_count,
		float **uvs, uint32_t *uvs_count,
		float width, float height, float depth,
		uint32_t width_segments, uint32_t height_segments, uint32_t depth_segments
	);
	extern int circle_buffer_geometry(
		uint8_t options,
		uint32_t **indices, uint32_t *indices_count,
		float **positions, uint32_t *positions_count,
		float **normals, uint32_t *normals_count,
		float **uvs, uint32_t *uvs_count,
		float radius, float theta_start, float theta_length,
		uint32_t segments
	);
	extern int cylinder_buffer_geometry(
		uint8_t options,
		uint32_t **indices, uint32_t *indices_count,
		float **positions, uint32_t *positions_count,
		float **normals, uint32_t *normals_count,
		float **uvs, uint32_t *uvs_count,
		float radius_top, float radius_bottom, float height,
		int open_ended, float theta_start, float theta_length,
		uint32_t radial_segments, uint32_t height_segments
	);

#ifdef __cplusplus
}
#endif

#endif // GEOM3D_H

#ifdef GEOM3D_IMPLEMENTATION

#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static float to_radians(float degrees) {
	return degrees * (M_PI / 180.0);
}

static int init_counts(uint8_t options,
	uint32_t *indices_count, uint32_t indices_count_val,
	uint32_t *positions_count, uint32_t positions_count_val,
	uint32_t *normals_count, uint32_t normals_count_val,
	uint32_t *uvs_count, uint32_t uvs_count_val) {
	if (options & GEOM3D_GENERATE_INDICES) {
		if (!indices_count) return 0;
		*indices_count = indices_count_val;
	}
	if (options & GEOM3D_GENERATE_POSITIONS) {
		if (!positions_count) return 0;
		*positions_count =  positions_count_val;
	}
	if (options & GEOM3D_GENERATE_NORMALS) {
		if (!normals_count) return 0;
		*normals_count =  normals_count_val;
	}
	if (options & GEOM3D_GENERATE_UVS) {
		if (!uvs_count)  return 0;
		*uvs_count = uvs_count_val;
	}

	return 1;
}

static int alloc_arrays(uint8_t options,
	uint32_t **indices, uint32_t indices_count,
	float **positions, uint32_t positions_count,
	float **normals, uint32_t normals_count,
	float **uvs, uint32_t uvs_count) {

	if (options & GEOM3D_GENERATE_INDICES) {
		if (!indices) return 0;
		*indices = (uint32_t*) calloc(indices_count, sizeof(uint32_t));
		if (!(*indices)) {
			return 0;
		}
	}
	if (options & GEOM3D_GENERATE_POSITIONS) {
		if (!positions) return 0;
		*positions = (float*) calloc(positions_count, sizeof(float));
		if (!(*positions)) {
			return 0;
		}
	}

	if (options & GEOM3D_GENERATE_NORMALS) {
		if (!normals) return 0;
		*normals = (float*) calloc(normals_count, sizeof(float));
		if (!(*normals)) {
			return 0;
		}
	}

	if (options & GEOM3D_GENERATE_UVS) {
		if (!uvs) return 0;
		*uvs = (float*) calloc(uvs_count, sizeof(float));
		if (!(*uvs)) {
			return 0;
		}
	}

	return 1;
}

int box_buffer_geometry(
	uint8_t options,
	uint32_t **indices, uint32_t *indices_count,
	float **positions, uint32_t *positions_count,
	float **normals, uint32_t *normals_count,
	float **uvs, uint32_t *uvs_count,
	float width, float height, float depth,
	uint32_t width_segments, uint32_t height_segments, uint32_t depth_segments
) {
	width_segments = width_segments < 1 ? 1 : width_segments;
	height_segments = height_segments < 1 ? 1 : height_segments;
	depth_segments = width_segments < 1 ? 1 : depth_segments;

	uint32_t points_count =
		(width_segments + 1) * (height_segments + 1) * 2 +
		(width_segments + 1) * (depth_segments  + 1) * 2 +
		(depth_segments + 1) * (height_segments + 1) * 2;
	uint32_t index_count = (
		width_segments * height_segments * 2 +
		width_segments * depth_segments * 2 +
		depth_segments * height_segments * 2) * 6;

	if (!init_counts(
		options,
		indices_count, index_count,
		positions_count, points_count * 3,
		normals_count, points_count * 3,
		uvs_count, points_count * 2)
	) return 0;

	if (!alloc_arrays(
		options,
		indices, *indices_count,
		positions, *positions_count,
		normals, *normals_count,
		uvs, *uvs_count)
	) return 0;

	struct panel { uint8_t u, v, w; int8_t udir, vdir; float width, height, depth; uint32_t gridx, gridy; };
	const struct panel panels[] = {
		{ 2, 1, 0, -1, -1, depth, height,  width,  depth_segments, height_segments },
		{ 2, 1, 0,  1, -1, depth, height, -width,  depth_segments, height_segments },
		{ 0, 2, 1,  1,  1, width, depth,   height, width_segments, depth_segments  },
		{ 0, 2, 1,  1, -1, width, depth,  -height, width_segments, depth_segments  },
		{ 0, 1, 2,  1, -1, width, height,  depth,  width_segments, height_segments },
		{ 0, 1, 2, -1, -1, width, height,  depth,  width_segments, height_segments }
	};

	uint32_t vertex_buffer_offset = 0, uv_buffer_offset = 0, index_buffer_offset = 0, number_of_vertices = 0;

	for (uint8_t i = 0; i < 6; i++) {
		float segment_width = panels[i].width / panels[i].gridx;
		float segment_height = panels[i].height / panels[i].gridy;

		float width_half = panels[i].width / 2;
		float height_half = panels[i].height / 2;
		float depth_half = panels[i].depth / 2;

		uint32_t gridx = panels[i].gridx + 1;
		uint32_t gridy = panels[i].gridy + 1;

		uint32_t vertex_counter = 0;

		float vector[3] = { 0, 0, 0 };

		for (uint32_t iy = 0; iy < gridy; iy++) {
			float y = 0;
			if (options & GEOM3D_GENERATE_POSITIONS)
				y = iy * segment_height - height_half;

			for (uint32_t ix = 0; ix < gridx; ix++) {
				float x = 0;
				if (options & GEOM3D_GENERATE_POSITIONS) {
					x = ix * segment_width - width_half;
					vector[panels[i].u] = x * panels[i].udir;
					vector[panels[i].v] = y * panels[i].vdir;
					vector[panels[i].w] = depth_half;

					(*positions)[vertex_buffer_offset]     = vector[0];
					(*positions)[vertex_buffer_offset + 1] = vector[1];
					(*positions)[vertex_buffer_offset + 2] = vector[2];
				}

				if (options & GEOM3D_GENERATE_NORMALS) {
					vector[panels[i].u] = 0;
					vector[panels[i].v] = 0;
					vector[panels[i].w] = panels[i].depth > 0 ? 1 : -1;

					(*normals)[vertex_buffer_offset]     = vector[0];
					(*normals)[vertex_buffer_offset + 1] = vector[1];
					(*normals)[vertex_buffer_offset + 2] = vector[2];
				}

				if (options & GEOM3D_GENERATE_UVS) {
					(*uvs)[uv_buffer_offset]     = ((float) ix) / gridx;
					(*uvs)[uv_buffer_offset + 1] = 1.0 - ((float) iy) / gridy;
				}

				vertex_buffer_offset += 3;
				uv_buffer_offset += 2;
				vertex_counter++;
			}
		}

		if (options & GEOM3D_GENERATE_INDICES) {
			for (uint32_t iy = 0; iy < panels[i].gridy; iy++) {
				for (uint32_t ix = 0; ix < panels[i].gridx; ix++) {
					uint32_t a = number_of_vertices + ix + gridx * iy;
					uint32_t b = number_of_vertices + ix + gridx * (iy + 1);
					uint32_t c = number_of_vertices + (ix + 1) + gridx * (iy + 1);
					uint32_t d = number_of_vertices + (ix + 1) + gridx * iy;

					(*indices)[ index_buffer_offset]     = a;
					(*indices)[ index_buffer_offset + 1] = b;
					(*indices)[ index_buffer_offset + 2] = d;

					(*indices)[index_buffer_offset + 3] = b;
					(*indices)[index_buffer_offset + 4] = c;
					(*indices)[index_buffer_offset + 5] = d;

					index_buffer_offset += 6;
				}
			}
			number_of_vertices += vertex_counter;
		}
	}

	return 1;
}

int circle_buffer_geometry(
	uint8_t options,
	uint32_t **indices, uint32_t *indices_count,
	float **positions, uint32_t *positions_count,
	float **normals, uint32_t *normals_count,
	float **uvs, uint32_t *uvs_count,
	float radius, float theta_start, float theta_length,
	uint32_t segments
) {
	segments = segments < 3 ? 3 : segments;
	uint32_t points_count = segments + 2;

	if (!init_counts(options, indices_count, 3 * segments, positions_count, points_count * 3, normals_count, points_count * 3, uvs_count, points_count * 2)) return 0;
	if (!alloc_arrays(options, indices, *indices_count, positions, *positions_count, normals, *normals_count, uvs, *uvs_count)) return 0;

	(*uvs)[0] = 0.5;
	(*uvs)[1] = 0.5;
	(*normals)[2] = 1.0;

	theta_start = to_radians(theta_start);
	theta_length = to_radians(theta_length);

	for (uint32_t s = 0, i = 3, j = 2; s <= segments; s++, i += 3, j += 2) {

		float pos1 = 0;
		float pos2 = 0;
		if ((options & GEOM3D_GENERATE_POSITIONS) || (options & GEOM3D_GENERATE_UVS)) {
			float segment = theta_start + ((float) s / segments) * theta_length;
			pos1 = radius * cos(segment);
			pos2 = radius * sin(segment);
		}

		if (options & GEOM3D_GENERATE_POSITIONS) {
			(*positions)[i] = pos1;
			(*positions)[i + 1] = pos2;
		}

		if (options & GEOM3D_GENERATE_NORMALS) {
			(*normals)[i + 2] = 1.0; // normal z
        }

		if (options & GEOM3D_GENERATE_UVS) {
			(*uvs)[j] = (pos1 / radius + 1.0) / 2.0;
			(*uvs)[j + 1] = (pos2 / radius + 1.0) / 2.0;
		}
	}

	if (options & GEOM3D_GENERATE_INDICES) {
		uint32_t *iter = *indices;
		for (uint32_t i = 1; i <= segments; i++) {
			*iter++ = i;
			*iter++ = i + 1;
			*iter++ = 0;
		}
	}

	return 1;
}

#endif // GEOM3D_IMPLEMENTATION
