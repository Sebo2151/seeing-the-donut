#pragma once

#include "util/Vectors.h"

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
};

