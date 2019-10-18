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

	void Draw();

private:
	// Computational data

	static const int grid_num_cols = 100;
	static const int grid_num_rows = 100;
	std::complex<float> tau;

	std::vector<Vector4> untransformed_verts;
	std::vector<int> indices;

	Vector4 translation;
	Matrix4 rotation;
	
	// GL data
	GLuint vert_shader_handle;
	GLuint frag_shader_handle;

	// TO DO: GL buffers


	void init_mesh();
	void init_shaders();

};

