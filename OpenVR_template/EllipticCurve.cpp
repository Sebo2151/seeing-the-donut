#include "EllipticCurve.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include "misc_util.h"

typedef  std::complex<double> C;

const double pi = 4.0 * atan(1.0);
const C ii = C(0.0, 1.0);

C theta(C z, C tau)
{
	C sum = 1;
	C nome = std::exp(pi * ii * tau);
	const int term_bound = 5;
	for (int n = 1; n <= term_bound; n++)
	{
		sum += 2.0 * pow(nome, n * n) * cos(2 * pi * n * z);
	}

	return sum;
}

C theta_01(C z, C tau)
{
	return theta(z + C(0.5, 0), tau);
}

C theta_10(C z, C tau)
{
	return exp(pi/4 * ii * tau + pi * ii * z) * theta(z + C(0.5, 0)*tau, tau);
}

C theta_11(C z, C tau)
{
	return exp(pi / 4 * ii * tau + pi * ii * (z + C(0.5, 0)))
		   * theta(z + C(0.5, 0) * tau + C(0.5, 0), tau);
}

C weierstrass_p(C z, C tau)
{
	C theta_0tau = theta(0, tau);
	C theta10_0tau = theta_10(0, tau);
	C theta01_ztau = theta_01(z, tau);
	C theta11_ztau = theta_11(z, tau);
	C term1 = pi * pi * theta_0tau * theta_0tau * theta10_0tau * theta10_0tau * theta01_ztau * theta01_ztau
		      / (theta11_ztau * theta11_ztau);
	C term2 = (pi * pi / 3.0) * (theta_0tau * theta_0tau * theta_0tau * theta_0tau
		+ theta10_0tau * theta10_0tau * theta10_0tau * theta10_0tau);

	return term1 - term2;
}

C weierstrass_p_derivative(C z, C tau)
{
	C sum = 0;

	const int sum_bounds = 10;
	for (int m = -sum_bounds; m <= sum_bounds; m++)
	{
		for (int n = -sum_bounds; n <= sum_bounds; n++)
		{
			C d = z + C(m,0) + (double)n * tau;
			sum += -2.0 / (d * d * d);
		}
	}

	return sum;
}

EllipticCurve::EllipticCurve()
{

	translation = Vector4(0, 0, 0, 0);
	rotation.identity();

	init_gl_info();
	init_mesh();
}

void EllipticCurve::Draw(Matrix4 projection)
{
	//dprintf("Here 2.5\n");
    
	// Build the vertex data to send to OpenGL
	std::vector<float> vert_data;
	std::vector<unsigned int> good_indices;

	for (int i = 0; i < untransformed_verts.size(); i++)
	{
		Vector4 transformed_vert = rotation * untransformed_verts[i] + translation;

		// Only copy 3 coords to project away a coordinate (namely, im(x), if no transformation done yet.)
		vert_data.push_back(transformed_vert.x/2);
		vert_data.push_back(transformed_vert.z/2);
		vert_data.push_back(transformed_vert.w/2);
		

		// Next is the color parameter, a # between 0 and 1 used for coloring.
		// Intended to make the projected-away coordinate somewhat visible.
		float color_param = atan(transformed_vert.y) / pi + 0.5f;
		vert_data.push_back(color_param );

		// Finally, the normal vector... search out some connected vertices and compute a cross product.
		// May be improved by computing more of these and taking an average.
		int this_m = i % grid_num_rows;
		int this_n = i / grid_num_rows;
		int next_m = (this_m + 1) % grid_num_cols;
		int next_n = (this_n + 1) % grid_num_rows;

		Vector4 v_up = rotation * untransformed_verts[this_m + next_n * grid_num_rows] + translation;
		Vector4 v_right = rotation * untransformed_verts[next_m + next_n * grid_num_rows] + translation;

		Vector3 v_up_3 = Vector3(v_up.x, v_up.z, v_up.w);
		Vector3 v_right_3 = Vector3(v_right.x, v_right.z, v_right.w);
		Vector3 v_3 = Vector3(transformed_vert.x, transformed_vert.z, transformed_vert.w);

		Vector3 normal = v_right_3 - v_3;
		normal = normal.cross(v_up_3 - v_3);
		normal.normalize();

		vert_data.push_back(normal.x);
		vert_data.push_back(normal.y);
		vert_data.push_back(normal.z);
	}

	const float cutoff_len = 30.0f;

	for (int m = 0; m < grid_num_cols; m++)
	{
		for (int n = 0; n < grid_num_rows; n++)
		{
			int i = 6 * (m * grid_num_rows + n);

			Vector3 tvs[6];
			for (int k = 0; k < 6; k++)
			{
				Vector4 tv = rotation * untransformed_verts[indices[i + k]] + translation;
				tvs[k] = Vector3(tv.x, tv.y, tv.z);
			}

			if ((tvs[0] - tvs[1]).length() < cutoff_len && (tvs[2] - tvs[1]).length() < cutoff_len && (tvs[2] - tvs[0]).length() < cutoff_len)
			{
				good_indices.push_back(indices[i]);
				good_indices.push_back(indices[i + 1]);
				good_indices.push_back(indices[i + 2]);
			}

			if ((tvs[3] - tvs[4]).length() < cutoff_len && (tvs[5] - tvs[4]).length() < cutoff_len && (tvs[5] - tvs[3]).length() < cutoff_len)
			{
				good_indices.push_back(indices[i + 3]);
				good_indices.push_back(indices[i + 4]);
				good_indices.push_back(indices[i + 5]);
			}
		}
	}

	
	// Now configure OpenGL to send stuff, then send the stuff.
	glUseProgram(gl_program_handle);
	glUniformMatrix4fv(gl_matrix_location, 1, GL_FALSE, projection.get());
	glBindVertexArray(gl_VAO_handle);

	//dprintf("Here 3\n");
	glBindBuffer(GL_ARRAY_BUFFER, gl_vert_data_buffer_handle);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float) * vert_data.size(), &vert_data[0], GL_DYNAMIC_DRAW);
	//dprintf("Here 4\n");

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gl_index_buffer_handle);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * good_indices.size(), &good_indices[0], GL_DYNAMIC_DRAW);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gl_index_buffer_handle);
	
	glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);

	glUseProgram(0);
}

void EllipticCurve::init_mesh()
{
	C tau = C(0, 1);

	// Compute the vertices!
	for (int i = 0; i < grid_num_cols; i++)
	{
		for (int j = 0; j < grid_num_rows; j++)
		{
			C z = ((i + 0.5)/grid_num_cols) * C(1, 0) + ((j + 0.5)/grid_num_rows) * tau;

			C p = weierstrass_p(z, tau);
			C p_prime = weierstrass_p_derivative(z, tau);

			untransformed_verts.push_back(Vector4(p.real(), p.imag(), p_prime.real(), p_prime.imag()));
		}
		dprintf("Progress on mesh: %i\n", i);
	}

	for (unsigned int i = 0; i < grid_num_cols; i++)
	{
		for (unsigned int j = 0; j < grid_num_rows; j++)
		{
			unsigned int next_i = (i + 1) % grid_num_cols;
			unsigned int next_j = (j + 1) % grid_num_rows;

			indices.push_back(j * grid_num_rows + i);
			indices.push_back(j * grid_num_rows + next_i);
			indices.push_back(next_j * grid_num_rows + next_i);

			indices.push_back(j * grid_num_rows + i);
			indices.push_back(next_j * grid_num_rows + next_i);
			indices.push_back(next_j * grid_num_rows + i);
		}
		
	}
}

void EllipticCurve::init_gl_info()
{
	gl_program_handle = CompileGLShader(
		"EllipticCurve",

		"#version 410\n"
		"uniform mat4 matrix;\n"
		"layout(location = 0) in vec3 position;\n"
		"layout(location = 1) in float color_param_in;\n"
		"layout(location = 2) in vec3 normal_in;\n"
		"out float color_param;\n"
		"out vec3 normal;\n"
		"void main()\n"
		"{\n"
		"    color_param = color_param_in;\n"
		"    normal = normal_in;\n"
		"    gl_Position = matrix * vec4(position, 1);\n"
		"}\n",

		"#version 410 core\n"
		"in float color_param;\n"
		"in vec3 normal;\n"
		"out vec3 output_color;\n"
		"void main()\n"
		"{\n"
		"   vec3 light1_dir = normalize(vec3(3,2,1));\n"
		"   vec3 light2_dir = normalize(vec3(-1, -1, 0));\n"
		"   vec3 diffuse = (color_param * vec3(247.0f/255.0f, 101.0f/255.0f, 15.0f/255.0f)\n"
		"                   + (1 - color_param) * vec3(31.0f/255.0f, 38.0f/255.0f, 244.0f/255.0f));\n"
		"   output_color = diffuse * ( max(0, dot(normal, light1_dir)) + max(0, 0.5f * dot(normal, light2_dir)));\n"
		"}\n"
	);

	if (gl_program_handle == 0)
	{
		dprintf("Fuck. Elliptic curve shader didn't compile.\n");
		exit(-1);
	}

	gl_matrix_location = glGetUniformLocation(gl_program_handle, "matrix");

	if (gl_matrix_location)
	{
		dprintf("The elliptic curve uniform matrix location is fucked.\n");
		exit(-1);
	}

	glGenVertexArrays(1, &gl_VAO_handle);
	glBindVertexArray(gl_VAO_handle);

	glGenBuffers(1, &gl_vert_data_buffer_handle);
	glBindBuffer(GL_ARRAY_BUFFER, gl_vert_data_buffer_handle);

	GLsizei stride = sizeof(float) * (3 + 1 + 3); // vertex + color param + normal

	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride, 0);

	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, stride, (const void*)(sizeof(float)*3));

	glEnableVertexAttribArray(2);
	glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, stride, (const void*)(sizeof(float)*4));

	//dprintf("Here 1\n");
	glGenBuffers(1, &gl_index_buffer_handle);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gl_index_buffer_handle);
	//dprintf("Here 2\n");

	glBindVertexArray(0);
}
