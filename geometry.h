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

//include implementation file
#include "geometry.tpp"