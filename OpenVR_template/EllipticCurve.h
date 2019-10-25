#pragma once

#include <complex>
#include <vector>

#include <GL/glew.h>
#include <SDL_opengl.h>
#include <GL/GLU.h>

#include "util/Vectors.h"
#include "util/Matrices.h"

class EllipticCurve
{
public:
	EllipticCurve();

	void Draw(Matrix4 projection);

private:
	// Computational data

	static const int grid_num_cols = 100;
	static const int grid_num_rows = 100;
	std::complex<float> tau;

	std::vector<Vector4> untransformed_verts;
	std::vector<unsigned int> indices;

	Vector4 translation;
	Matrix4 rotation;
	
	// GL data
	GLuint gl_program_handle = 0;
	GLuint gl_matrix_location = 0;

	GLuint gl_vert_data_buffer_handle = 0;
	GLuint gl_index_buffer_handle = 0;
	GLuint gl_VAO_handle = 0;

	// TO DO: GL buffers

	void init_mesh();
	void init_gl_info();

};

