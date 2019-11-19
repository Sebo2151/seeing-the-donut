#pragma once

#include "util/Vectors.h"
#include "util/Matrices.h"

class Quat
{
public:
	float x, y, z, w; // xi + yj + zk + w

	Quat(Vector4 in)
	{
		x = in.x;
		y = in.y;
		z = in.z;
		w = in.w;
	}

	Quat(float x, float y, float z, float w)
	{
		this->x = x;
		this->y = y;
		this->z = z;
		this->w = w;
	}

	float norm()
	{
		return sqrt(x * x + y * y + z * z + w * w);
	}

	Quat conj()
	{
		return Quat(-x, -y, -z, w);
	}

	Quat inv()
	{
		float n_sqr = x * x + y * y + z * z + w * w;
		return Quat(-x / n_sqr, -y / n_sqr, -z / n_sqr, w / n_sqr);
	}

	Quat normalized()
	{
		float n = this->norm();
		return Quat(x / n, y / n, z / n, w / n);
	}

	Vector4 toVector4()
	{
		return Vector4(x, y, z, w);
	}

	Quat& operator+(const Quat& other) const
	{
		return Quat(x + other.x, y + other.y, z + other.z, w + other.w);
	}

	Quat& operator-(const Quat& other) const
	{
		return Quat(x - other.x, y - other.y, z - other.z, w - other.w);
	}

	Quat& operator*(const Quat& r) const
	{
		return Quat(
			x*r.w + y*r.z - z*r.y + w*r.x,
			-x*r.z + y*r.w + z*r.x + w*r.y,
			x*r.y - y*r.x + z*r.w + w*r.z,
			w* r.w - x * r.x - y * r.y - z * r.z
		);
	}

	Quat& operator/(float d)
	{
		return Quat(x / d, y / d, z / d, w / d);
	}

	static Quat fromMat3(Matrix3 mat)
	{
		// Adjusted from http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
		const float* m = mat.get();

		float tr = m[0*3 + 0] + m[1*3 + 1] + m[2*3 + 2];
		float w, x, y, z;

		if (tr > 0) {
			float S = sqrt(tr + 1.0) * 2; // S=4*qw 
			w = 0.25 * S;
			x = (m[2*3 + 1] - m[1*3 + 2]) / S;
			y = (m[0*3 + 2] - m[2*3 + 0]) / S;
			z = (m[1*3 + 0] - m[0*3 + 1]) / S;
		}
		else if ((m[0*3 + 0] > m[1*3 + 1]) && (m[0*3 + 0] > m[2*3 + 2])) {
			float S = sqrt(1.0 + m[0*3 + 0] - m[1*3 + 1] - m[2*3 + 2]) * 2; // S=4*qx 
			w = (m[2*3 + 1] - m[1*3 + 2]) / S;
			x = 0.25 * S;
			y = (m[0*3 + 1] + m[1*3 + 0]) / S;
			z = (m[0*3 + 2] + m[2*3 + 0]) / S;
		}
		else if (m[1*3 + 1] > m[2*3 + 2]) {
			float S = sqrt(1.0 + m[1*3 + 1] - m[0*3 + 0] - m[2*3 + 2]) * 2; // S=4*qy
		    w = (m[0*3 + 2] - m[2*3 + 0]) / S;
			x = (m[0*3 + 1] + m[1*3 + 0]) / S;
			y = 0.25 * S;
			z = (m[1*3 + 2] + m[2*3 + 1]) / S;
		}
		else {
			float S = sqrt(1.0 + m[2*3 + 2] - m[0*3 + 0] - m[1*3 + 1]) * 2; // S=4*qz
			w = (m[1*3 + 0] - m[0*3 + 1]) / S;
			x = (m[0*3 + 2] + m[2*3 + 0]) / S;
			y = (m[1*3 + 2] + m[2*3 + 1]) / S;
			z = 0.25 * S;
		}

		return Quat(x, y, z, w);
	}

	static Matrix4 toRot(Quat first, Quat second)
	{
		// Might not be efficient... but it should work!
		Quat second_inv = second.inv();
		Quat col1 = first * Quat(1, 0, 0, 0) * second_inv;
		Quat col2 = first * Quat(0, 1, 0, 0) * second_inv;
		Quat col3 = first * Quat(0, 0, 1, 0) * second_inv;
		Quat col4 = first * Quat(0, 0, 0, 1) * second_inv;

		return Matrix4(
			col1.x, col2.x, col3.x, col4.x,
			col1.y, col2.y, col3.y, col4.y,
			col1.z, col2.z, col3.z, col4.z,
			col1.w, col2.w, col3.w, col4.w);
	}
};

