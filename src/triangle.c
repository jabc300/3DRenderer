#include "display.h"
#include "swap.h"
#include "triangle.h"

vec3_t get_triangle_normal(vec4_t vertices[3]) {
	// Check backface culling
	vec3_t vector_a = vec3_from_vec4(vertices[0]);
	vec3_t vector_b = vec3_from_vec4(vertices[1]);
	vec3_t vector_c = vec3_from_vec4(vertices[2]);
			
	// Get vector from A to B and A to C
	vec3_t vector_ab = vec3_sub(vector_b, vector_a); 
	vec3_t vector_ac = vec3_sub(vector_c, vector_a); 
	vec3_normalize(&vector_ab);
	vec3_normalize(&vector_ac);

	// Compute the face normal (using cross product to find perpendicular)
	vec3_t normal = vec3_cross(vector_ab, vector_ac);
			
	// Normalize the face normal vector
	vec3_normalize(&normal);

	return normal;
}

void draw_triangle_pixel(int x, int y, uint32_t color, vec4_t point_a, vec4_t point_b, vec4_t point_c) {

	vec2_t p = { x, y };
	vec2_t a = vec2_from_vec4(point_a);
	vec2_t b = vec2_from_vec4(point_b);
	vec2_t c = vec2_from_vec4(point_c);

	vec3_t weights = barycentric_weights(a, b, c, p);
	
	float alpha = weights.x;
	float beta = weights.y;
	float gamma = weights.z;

	// Interpolate the value of 1/w for the current pixel
	float interpolated_reciprocal_w = (1 / point_a.w) * alpha + (1 / point_b.w) * beta + (1 / point_c.w) * gamma;

	// Adjust 1/w so pixels closer to the camera have smaller values. 
	interpolated_reciprocal_w = 1.0 - interpolated_reciprocal_w;

	// Only draw pixel if the depth value is less than the one previously store in the z buffer
	if(interpolated_reciprocal_w < get_zbuffer_at(x, y)) {	
		draw_pixel(x, y, color);

		update_zbuffer_at(x, y, interpolated_reciprocal_w);
	}
}

void draw_filled_triangle(
	int x0, int y0, float z0, float w0, 
	int x1, int y1, float z1, float w1,
	int x2, int y2, float z2, float w2, 
	uint32_t color) {
	// We need to sort the vertices by y-coordinate ascending (y0 < y1 < y2)
	if (y0 > y1) {
		int_swap(&y0, &y1);
		int_swap(&x0, &x1);
		float_swap(&z0, &z1);
		float_swap(&w0, &w1);
	}
	
	if (y1 > y2) {
		int_swap(&y1, &y2);
		int_swap(&x1, &x2);
		float_swap(&z1, &z2);
		float_swap(&w1, &w2);
	}

	if (y0 > y1) {
		int_swap(&y0, &y1);
		int_swap(&x0, &x1);
		float_swap(&z0, &z1);
		float_swap(&w0, &w1);
	}
	
	// Create vector points and texture coordinates after we sort the vertices
	vec4_t point_a = { x0, y0, z0, w0 };
	vec4_t point_b = { x1, y1, z1, w1 };
	vec4_t point_c = { x2, y2, z2, w2 };

	// Render the upper part of the triangle (flat-bottom)
	float inv_slope1 = 0;
	float inv_slope2 = 0;

	if (y1 - y0 != 0) inv_slope1 = (float)(x1 - x0) / abs(y1 - y0);
	if (y2 - y0 != 0) inv_slope2 = (float)(x2 - x0) / abs(y2 - y0);

	if(y1 - y0 != 0) {
		for (int y = y0; y <= y1; y++) {
			int x_start = x1 + (y - y1) * inv_slope1;
			int x_end = x0 + (y - y0) * inv_slope2;
				
			if (x_end < x_start) {
				int_swap(&x_start, &x_end);
			}

			for (int x = x_start; x < x_end; x++) 
			{
				draw_triangle_pixel(x, y, color, point_a, point_b, point_c);
			}

		}
	}
	// Render the bottom part of the triangle (flat-top)
	inv_slope1 = 0;
	inv_slope2 = 0;


	if (y2 - y1 != 0) inv_slope1 = (float)(x2 - x1) / abs(y2 - y1);
	if (y2 - y0 != 0) inv_slope2 = (float)(x2 - x0) / abs(y2 - y0);
	
	if (y2 - y1 != 0) {
		for (int y = y1; y <= y2; y++) {
			int x_start = x1 + (y - y1) * inv_slope1;
			int x_end = x0 + (y - y0) * inv_slope2;
		
			if (x_end < x_start) {
				int_swap(&x_start, &x_end);
			}

			for (int x = x_start; x < x_end; x++) {
				draw_triangle_pixel(x, y, color, point_a, point_b, point_c);
			}
		}
	}
	
}

void draw_triangle(int x0, int y0, int x1, int y1, int x2, int y2, uint32_t color) {
	draw_line(x0, y0, x1, y1, color);
	draw_line(x1, y1, x2, y2, color);
	draw_line(x2, y2, x0, y0, color);
}

vec3_t barycentric_weights(vec2_t a, vec2_t b, vec2_t c, vec2_t p) {
	vec2_t ac = vec2_sub(c, a);
	vec2_t ab = vec2_sub(b, a);
	vec2_t ap = vec2_sub(p, a);
	vec2_t pc = vec2_sub(c, p);
	vec2_t pb = vec2_sub(b, p);

	float area_parallelogram_abc = (ac.x * ab.y - ac.y * ab.x);

	float alpha = (pc.x * pb.y - pc.y * pb.x) / area_parallelogram_abc;
	
	float beta = (ac.x * ap.y - ac.y * ap.x) / area_parallelogram_abc;

	float gamma = 1 - alpha - beta;

	vec3_t weights = { alpha, beta, gamma };
	return weights;
}

void draw_texel(
	int x, int y, upng_t* texture,
	vec4_t point_a, vec4_t point_b, vec4_t point_c,
	tex2_t a_uv, tex2_t b_uv, tex2_t c_uv
) {	
	vec2_t p = { x, y };
	vec2_t a = vec2_from_vec4(point_a);
	vec2_t b = vec2_from_vec4(point_b);
	vec2_t c = vec2_from_vec4(point_c);

	vec3_t weights = barycentric_weights(a, b, c, p);
	
	float alpha = weights.x;
	float beta = weights.y;
	float gamma = weights.z;

	// Perform interpolation of all u and v values using barycentric weights
	float interpolated_u;
	float interpolated_v;
	float interpolated_reciprocal_w;

	interpolated_u = (a_uv.u / point_a.w) * alpha + (b_uv.u / point_b.w) * beta + (c_uv.u / point_c.w) * gamma;
	interpolated_v = (a_uv.v / point_a.w) * alpha + (b_uv.v / point_b.w) * beta + (c_uv.v / point_c.w) * gamma;
	
	// Interpolate the value of 1/w for the current pixel
	interpolated_reciprocal_w = (1 / point_a.w) * alpha + (1 / point_b.w) * beta + (1 / point_c.w) * gamma;

	interpolated_u /= interpolated_reciprocal_w;
	interpolated_v /= interpolated_reciprocal_w;

	// Get the mesh texture width and height dimensions
	int texture_width = upng_get_width(texture);
	int texture_height = upng_get_height(texture);

	// Map the UV coordinates to the full texture width and height
	int tex_x = abs((int)(interpolated_u * texture_width));
	int tex_y = abs((int)(interpolated_v * texture_height));
	
	// Adjust 1/w so pixels closer to the camera have smaller values. 
	interpolated_reciprocal_w = 1.0 - interpolated_reciprocal_w;

	// Only draw pixel if the depth value is less than the one previously store in the z buffer
	if(interpolated_reciprocal_w < get_zbuffer_at(x, y)) {	
		// Get the buffer of colors from the texture
		uint32_t* texture_buffer = (uint32_t*)upng_get_buffer(texture);

		int texIndex = ((texture_width * tex_y) + tex_x) % (texture_width * texture_height); 
		draw_pixel(x, y, texture_buffer[texIndex]);

		update_zbuffer_at(x, y, interpolated_reciprocal_w);
	}
}

void draw_textured_triangle(
	int x0, int y0, float z0, float w0, float u0, float v0, 
	int x1, int y1, float z1, float w1, float u1, float v1, 
	int x2, int y2, float z2, float w2, float u2, float v2,
	upng_t* texture
) {
	if (y0 > y1) {
		int_swap(&y0, &y1);
		int_swap(&x0, &x1);
		float_swap(&z0, &z1);
		float_swap(&w0, &w1);
		float_swap(&u0, &u1);
		float_swap(&v0, &v1);
	}
	
	if (y1 > y2) {
		int_swap(&y1, &y2);
		int_swap(&x1, &x2);
		float_swap(&z1, &z2);
		float_swap(&w1, &w2);
		float_swap(&u1, &u2);
		float_swap(&v1, &v2);
	}

	if (y0 > y1) {
		int_swap(&y0, &y1);
		int_swap(&x0, &x1);
		float_swap(&z0, &z1);
		float_swap(&w0, &w1);
		float_swap(&u0, &u1);
		float_swap(&v0, &v1);
	}
	
	// Flip V component
	v0 = 1.0 - v0;
	v1 = 1.0 - v1;
	v2 = 1.0 - v2;

	// Create vector points and texture coordinates after we sort the vertices
	vec4_t point_a = { x0, y0, z0, w0 };
	vec4_t point_b = { x1, y1, z1, w1 };
	vec4_t point_c = { x2, y2, z2, w2 };
	tex2_t a_uv = { u0, v0 };
	tex2_t b_uv = { u1, v1 };
	tex2_t c_uv = { u2, v2 };

	// Render the upper part of the triangle (flat-bottom)
	float inv_slope1 = 0;
	float inv_slope2 = 0;

	if (y1 - y0 != 0) inv_slope1 = (float)(x1 - x0) / abs(y1 - y0);
	if (y2 - y0 != 0) inv_slope2 = (float)(x2 - x0) / abs(y2 - y0);

	if(y1 - y0 != 0) {
		for (int y = y0; y <= y1; y++) {
			int x_start = x1 + (y - y1) * inv_slope1;
			int x_end = x0 + (y - y0) * inv_slope2;
				
			if (x_end < x_start) {
				int_swap(&x_start, &x_end);
			}

			for (int x = x_start; x < x_end; x++) 
			{
				draw_texel(x, y, texture, point_a, point_b, point_c, a_uv, b_uv, c_uv);
			}

		}
	}
	// Render the bottom part of the triangle (flat-top)
	inv_slope1 = 0;
	inv_slope2 = 0;


	if (y2 - y1 != 0) inv_slope1 = (float)(x2 - x1) / abs(y2 - y1);
	if (y2 - y0 != 0) inv_slope2 = (float)(x2 - x0) / abs(y2 - y0);
	
	if (y2 - y1 != 0) {
		for (int y = y1; y <= y2; y++) {
			int x_start = x1 + (y - y1) * inv_slope1;
			int x_end = x0 + (y - y0) * inv_slope2;
		
			if (x_end < x_start) {
				int_swap(&x_start, &x_end);
			}

			for (int x = x_start; x < x_end; x++) {
				draw_texel(x, y, texture, point_a, point_b, point_c, a_uv, b_uv, c_uv);
			}
		}
	}
}
