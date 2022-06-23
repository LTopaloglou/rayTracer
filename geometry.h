#include <cstdlib>
#include <cmath>
//This file contains all basic geometric classes necessary for ray tracing


template <typename T> class Vector3 {
    //Vector3 is a vector in three dimensions
public:
    //--Class variables--
    //x y and z components
    T x, y, z;

    //--Constructors--
    //default constructor setting components to 0
    Vector3();
    //parametric constructor
    Vector3(T x_, T y_, T z_);

    //--Accessors--
    //overloaded index operator so vector can be accessed like an array
    //this one lets components be seen but not changed
    T operator[](int i) const;
    //this one lets components be changed because it returns reference to them
    T &operator[](int i);
    
    //--Vector operations--
    //vector adition and +=.
    Vector3<T> operator+ (const Vector3<T> &rhs) const;
    Vector3<T>& operator+= (const Vector3<T> &rhs);

    //Vector subtraction & negative vector
    Vector3<T> operator-() const;
    Vector3<T> operator-(const Vector3<T> &rhs) const;
    Vector3<T>& operator-=(const Vector3<T> &rhs);

    //scalar multiplication and *=
    Vector3<T> operator* (T i) const;
    Vector3<T>& operator*= (T i);

    //scalar division and /=
    Vector3<T> operator/ (T i) const;
    Vector3<T>& operator/= (T i);

    //dot product defined as the * operator
    T operator* (const Vector3<T> &rhs) const;

    //absolute value of a vector
    Vector3<T> abs() const;

    //magnitude of a vector
    float Mag() const;

};

//Scalar multiplication in the case of datatype T * Vector3
template <typename T>
Vector3<T> operator* (T i, const Vector3<T> &vec);

//Dot product defined as Dot(vec1, vec2)
template <typename T>
T Dot(const Vector3<T> &lhs, const Vector3<T> &rhs);

//Cross product
template <typename T>
Vector3<T> Cross(const Vector3<T> &lhs, const Vector3<T> &rhs);

//Unit vector aka normalized vector
template <typename T>
Vector3<T> Unit(const Vector3<T> &vec);

//shortened name for vectors & zero vector
typedef Vector3<int> Vector3i;
typedef Vector3<float> Vector3f;
const Vector3<int> zero3i();
const Vector3<float> zero3f();

template <typename T> class Point3 {
    //Point3 is a point in 3 dimensions
public:
    //--Class variables--
    //x, y, z coordinates
    T x, y, z;

    //--Constructors--
    //Default constructor
    Point3();
    //Parametric constructor
    Point3(T x_, T y_, T z_);

    //--Accessors--
    //overloaded index operator to allow point to be accesssed like an array
    T operator[] (int i) const;
    T& operator[](int i);

    //--Point operations--
    
    //Point-vector addition and +=
    //Note: + and - only work when its Point + Vector (this returns a point)
    //Vector + point would not make sense since it would imply a vector is returned
    //and points cannot be added to vectors
    Point3<T> operator+ (const Vector3<T> &vec) const;
    Point3<T>& operator+= (const Vector3<T> &vec);
    //Point-vector subtraction
    Point3<T> operator- (const Vector3<T> &vec) const;
    Point3<T>& operator-= (const Vector3<T> &vec);

    //Point-point subtraction
    //e.g. Vector AB = point B - point A
    Vector3<T> operator- (const Point3<T> &rhs) const;

};

//Distance between two points
template <typename T>
float Distance(const Point3<T> &lhs,const Point3<T> &rhs);

//Linear Interpolation between two points
template <typename T>
Point3<float> Lerp(const Point3<T> &lhs, const Point3<T> &rhs, float n);

//definitions to make life easy
typedef Point3<int> Point3i;
typedef Point3<float> Point3f;
template <typename T> class Normal3 {
    //class that represents a normal to a surface.
    //note that although similar to vector it is not the same
public:
    //--Class variables--
    T x, y, z;

    //--Constructors--
    //defining a normal from a vector
    Normal3(const Vector3<T> & vec);
};

class Ray {
    //this class represents a ray being shot from the camera
public:
    //--Class Variables--
    Point3f origin; // origin of the vector
    Vector3f direction; //direction vector of the ray
    mutable float max; //max scalar multiple of direction that is allowed
    float time; //time associated with the ray
    const Medium *medium; //the medium (e.g. air, fog, water etc) that the ray travels thrhough

    //--Constructors--
    //default constructor
    Ray();
    //origin - direction constructor
    Ray(const Point3f & origin_, const Vector3f & direction_);

    //Overloaded Operators
    Point3f operator()(float f);
};

template <typename T> class Bounds3 {
    //class that represents a bounding box of an object
public:
    //--Class variables--
    Point3<T> min, max;

    //--Constructors--
    //default constructor
    Bounds3();
    //min point + x,y,z, distance constructor
    Bounds3(const Point3<T> & min_, T x_, T y_, T z_);
    //two-point constructor. Note that since p2 is not always greater than p1, 
    //we must first work out the maximums in each x,y,z directions and create new points
    Bounds3(const Point3<T> & p1, const Point3<T> & p2);
    //single point bound (bounds nothing but that point)
    Bounds3(const Point3<T> & point);
};

class Matrix4x4 {
    //class that represents a 4x4 matrix. It will be used for both
    //linear transformations & transformations between coordinate systems
    //which is why it is 4x4 and not 3x3
public:
    //--class variables--
    float m[4][4];

    //--Constructors--
    //default constructor creates an identity matrix
    Matrix4x4();
    //paramatric constructor
    Matrix4x4(float m00, float m01, float m02, float m03,
              float m10, float m11, float m12, float m13,
              float m20, float m21, float m22, float m23,
              float m30, float m31, float m32, float m33);


    //Overloaded operators
    //matrix addition
    Matrix4x4 operator+(const Matrix4x4 & rhs);
};

class Medium {

};

//include implementation file
#include "geometry.tpp"