#include "EllipticCurve.h"

#include <cmath>
#include "misc_util.h"

typedef  std::complex<double> C;

const double pi = 2.0 * atan(1.0);
// No index i's in this file!!!
const C ii = C(0.0, 1.0);

C theta(C z, C tau)
{
	C sum = 0;
    
	const int term_bound = 25;
	for (int n = -term_bound; n <= term_bound; n++)
	{
		sum += std::exp(pi * n * n * ii * tau + 2 * pi * n * ii * z);
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
		   * theta(z + C(0.5, 0) * tau, tau);
}

C weierstrass_p(C z, C tau)
{
	C theta_0tau = theta(0, tau);
	C theta10_0tau = theta_10(0, tau);
	C theta01_ztau = theta_01(z, tau);
	C theta11_ztau = theta_11(z, tau);
	C term1 = pi * pi * theta10_0tau * theta10_0tau * theta01_ztau * theta01_ztau / (theta11_ztau * theta11_ztau);
	C term2 = (pi * pi / 3.0) * (theta_0tau * theta_0tau * theta_0tau * theta_0tau
		+ theta10_0tau * theta10_0tau * theta10_0tau * theta10_0tau);

	return term1 - term2;
}

C weierstrass_p_derivative(C z, C tau)
{
	C sum = 0;

	const int sum_bounds = 50;
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

	init_mesh();
	init_gl_info();
}

void EllipticCurve::Draw()
{
}

void EllipticCurve::init_mesh()
{
	C tau = C(0, 1);

	// Compute the vertices!
	for (int i = 0; i < grid_num_cols; i++)
	{
		for (int j = 0; j < grid_num_rows; j++)
		{
			C z = (i + 0.5) * C(1, 0) + (j + 0.5) * tau;

			C p = weierstrass_p(z, tau);
			C p_prime = weierstrass_p_derivative(z, tau);

			untransformed_verts.push_back(Vector4(p.real(), p.imag(), p_prime.real(), p_prime.imag()));
		}
	}

	for (int i = 0; i < grid_num_cols; i++)
	{
		for (int j = 0; j < grid_num_rows; j++)
		{
			int next_i = (i + 1) % grid_num_cols;
			int next_j = (j + 1) % grid_num_rows;

			indices.push_back(j * grid_num_cols + i);
			indices.push_back(j * grid_num_cols + next_i);
			indices.push_back(next_j * grid_num_cols + next_i);

			indices.push_back(j * grid_num_cols + i);
			indices.push_back(next_j * grid_num_cols + next_i);
			indices.push_back(next_j * grid_num_cols + i);
		}
	}
}

void EllipticCurve::init_gl_info()
{
	gl_program_handle = CompileGLShader(
		"EllipticCurve",

		"#version 410\n"
		"uniform mat4 matrix;\n"
		"layout(location = 0) in vec4 position;\n"
		"layout(location = 1) in float color_param_in;\n"
		"layout(location = 2) in vec3 normal_in;\n"
		"out float color_param;\n"
		"out vec3 normal;\n"
		"void main()\n"
		"{\n"
		"    color_param = color_param_in;\n"
		"    normal = normal_in;\n"
		"    gl_Position = matrix * position;\n"
		"}\n",

		"#version 410 core\n"
		"in float color_param;\n"
		"in vec3 normal;\n"
		"out vec3 output_color;\n"
		"void main()\n"
		"{\n"
		"   output_color = (color_param * vec3(83.0f/255.0f, 236.0f/255.0f, 210.0f/255.0f)"
		"                   + (1 - color_param) * vec3(49.0f/255.0f, 124.0f/255.0f, 238.0f/255.0f))"
		"                * (abs(dot(normal, vec3(.577,.577,.577))));\n"
		"}\n"
	);

	if (gl_program_handle == 0)
	{
		dprintf("Fuck. Elliptic curve shader didn't compile.\n");
	}

	gl_matrix_location = glGetUniformLocation(gl_program_handle, "matrix");

	if (gl_matrix_location)
	{
		dprintf("The elliptic curve uniform matrix location is fucked.\n");
	}


}
