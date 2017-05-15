// Transform.cpp: implementation of the Transform class.

// Note: when you construct a matrix using mat4() or mat3(), it will be COLUMN-MAJOR
// Keep this in mind in readfile.cpp and display.cpp
// See FAQ for more details or if you're having problems.

#include "Transform.h"

// Helper rotation function.  Please implement this.  
mat3 Transform::rotate(const float degrees, const vec3& axis) 
{
	//Shorthands for commonly used parameters
	float cos = cosf(degrees * pi / 180);
	float sin = sinf(degrees * pi / 180);
	vec3 normal = normalize(axis);
	float a = normal[0];
	float b = normal[1];
	float c = normal[2];

	//Axis-Angle formula matrices
	mat3 A(a*a, a*b, a*c, a*b, b*b, b*c, a*c, b*c, c*c);
	mat3 B(0, c, -b, -c, 0, a, b, -a, 0);

	//Apply the formula
	return mat3()*cos + A*(1 - cos) + B*sin;
}

void Transform::left(float degrees, vec3& eye, vec3& up) 
{
	eye = rotate(degrees, up)*eye;
}

void Transform::up(float degrees, vec3& eye, vec3& up) 
{
	//Find the orthogonal axis
	vec3 A = cross(eye, up);
	eye = rotate(degrees, A)*eye;
	up = rotate(degrees, A)*up;
}

mat4 Transform::lookAt(const vec3 &eye, const vec3 &center, const vec3 &up) 
{
	//Coordinate vectors
	vec3 w = normalize(eye);
	vec3 u = normalize(cross(up, w));
	vec3 v = normalize(cross(w, u));

	//Construct matrix
	mat4 A(u.x, v.x, w.x, 0, u.y, v.y, w.y, 0, u.z, v.z, w.z, 0, 0, 0, 0, 1);
	mat4 B(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, -eye.x, -eye.y, -eye.z, 1);
	return A*B;
}

mat4 Transform::perspective(float fovy, float aspect, float zNear, float zFar)
{
	//Get d,A,B from fovy,aspect,f,n
	float d = 1/tanf(fovy / 2 * pi / 180);
	float A = -(zFar + zNear) / (zFar - zNear);
	float B = -(2 * zFar*zNear) / (zFar - zNear);
    return mat4(d/aspect, 0, 0, 0,
				0, d, 0, 0,
				0, 0, A, -1,
				0, 0, B, 0);
}

mat4 Transform::scale(const float &sx, const float &sy, const float &sz) 
{
	return mat4(sx, 0, 0, 0,
				0, sy, 0, 0,
				0, 0, sz, 0,
				0, 0, 0, 1);
}

mat4 Transform::translate(const float &tx, const float &ty, const float &tz) 
{
	return mat4(1, 0, 0, 0,
				0, 1, 0, 0,
				0, 0, 1, 0,
				tx, ty, tz, 1);
}

// To normalize the up direction and construct a coordinate frame.  
// As discussed in the lecture.  May be relevant to create a properly 
// orthogonal and normalized up. 
// This function is provided as a helper, in case you want to use it. 
// Using this function (in readfile.cpp or display.cpp) is optional.  

vec3 Transform::upvector(const vec3 &up, const vec3 & zvec) 
{
    vec3 x = glm::cross(up,zvec); 
    vec3 y = glm::cross(zvec,x); 
    vec3 ret = glm::normalize(y); 
    return ret; 
}


Transform::Transform()
{

}

Transform::~Transform()
{

}