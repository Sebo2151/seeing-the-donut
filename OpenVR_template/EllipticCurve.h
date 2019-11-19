#pragma once

#include <complex>
#include <vector>

#include <GL/glew.h>
#include <SDL_opengl.h>
#include <GL/GLU.h>

#include <openvr.h>

#include "util/Vectors.h"
#include "util/Matrices.h"

class EllipticCurve
{
public:
	EllipticCurve();

	void Draw(Matrix4 projection);
	void Update(vr::InputDigitalActionData_t leftButtonActionData, vr::InputDigitalActionData_t rightButtonActionData,
		vr::InputPoseActionData_t leftPose, vr::InputPoseActionData_t rightPose);

private:
	// Computational data

	static const int grid_num_cols = 200;
	static const int grid_num_rows = 200;
	std::complex<float> tau;

	std::vector<Vector4> untransformed_verts;
	std::vector<unsigned int> indices;
	
	// GL data
	GLuint gl_program_handle = 0;
	GLuint gl_matrix_location = 0;

	GLuint gl_vert_data_buffer_handle = 0;
	GLuint gl_index_buffer_handle = 0;
	GLuint gl_VAO_handle = 0;

	// TO DO: GL buffers

	void init_mesh();
	void init_gl_info();

	Vector4 translation;
	Matrix4 rotation;

	// UI data
	Matrix3 storedRotLeft;
	Matrix3 storedRotRight;
	Matrix4 oldRotation;
	bool rotating;
};

