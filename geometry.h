//This file contains all basic geometric classes necessary for ray tracing


template <typename T> class Vector3 {
    //Vector3 is a vector in three dimensions
public:
    //x y and z components
    T x, y, z;
    //default constructor setting components to 0
    Vector3();
    //overloaded index operator so vector can be accessed like an array
    //this one lets components be seen but not changed
    T operator[](int i) const;
    //this one lets components be changed because it returns reference to them
    T &operator[](int i);

};

//shortened name for vectors & zero vector
typedef Vector3<int> Vector3i;
typedef Vector3<float> Vector3f;
const Vector3<int> zero3i();
const Vector3<float> zero3f();

//include implementation file
#include "geometry.tpp"