//#include "geometry.h"
//This is an implementation file because template classes cannot be implemented 
//in the .cpp and i want to keep the .h file clean

//Vector3 implementation
template <typename T> 
Vector3<T>::Vector3() {
    x = y = z = 0;
}

template <typename T> 
Vector3<T>::Vector3(T x_, T y_, T z_) {
    x = x_;
    y = y_;
    z = z_;
}

template <typename T>
T Vector3<T>::operator[] (int i) const {
    //check for invalid i, throw error if true
    Assert(i >= 0 && i <= 2);
    if (i == 0): return x;
    if (i == 1): return y;
    return z;
}

template <typename T>
T& Vector3<T>::operator[] (int i) {
    //check for invalid i, throw error if true
    Assert(i >= 0 && i <= 2);
    if (i == 0): return x;
    if (i == 1): return y;
    return z;
}

template <typename T>
Vector3<T> Vector3<T>::operator+ (const Vector3<T> &rhs) const {
    return Vector3(x + rhs.x, y + rhs.y, z + rhs.z);
}

template <typename T>
Vector3<T>& Vector3<T>::operator+= (const Vector3<T> &rhs) {
    x += rhs.x;
    y += rhs.y;
    z += rhs.z;
    return *this;
}

template <typename T>
Vector3<T> Vector3<T>::operator-() const{
    //this function returns the negative of the vector
    return Vector3(-x, -y, -z);
}

template <typename T>
Vector3<T> Vector3<T>::operator-(const Vector3<T> &rhs) const {
    return (x-rhs.x, y-rhs.y, z-rhs.z);
}

template <typename T>
Vector3<T>& Vector3<T>::operator-=(const Vector3<T> &rhs) {
    x -= rhs.x;
    y -= rhs.y;
    z -= rhs.z;
    return *this;
}

template <typename T>
Vector3<T> Vector3<T>::operator* (T i) const {
    return Vector(x*i, y*i, z*i);
}

template <typename T>
Vector3<T>& Vector3<T>::operator*= (T i) {
    x *= i;
    y *= i;
    z *= i;
    return *this;
}

template <typename T>
Vector3<T> operator* (T i, Vector3<T> vec) {
    return Vector(vec.x*i, vex.y*i, vec.z*i);
}

template <typename T>
Vector3<T>  Vector3<T>::operator/ (T i)const {
    //make sure were not dividing by zero
    Assert(i != 0);
    //calculate reciprocal of i first since this means only 1 division compared to 3
    //because division takes longer than multiplication
    float reciprocal = 1/i;
    return Vector(x*reciprocal, y*reciprocal, z*reciprocal);
}

template <typename T>
Vector3<T>& Vector3<T>::operator/= (T i) {
    //make sure were not dividing by zero
    Assert(i != 0);
    //calculate reciprocal of i first since this means only 1 division compared to 3
    //because division takes longer than multiplication
    float reciprocal = 1/i;
    x *= reciprocal;
    y *= reciprocal;
    z *= reciprocal;
    return *this;
}

template <typename T>
T Vector3<T>::operator* (const Vector3<T> &rhs) const {
    return x*rhs.x + y*rhs.y + z*rhs.z;
}

template <typename T>
T Dot(const Vector3<T> &lhs, const Vector3<T> &rhs) {
    return lhs.x*rhs.x + lhs.y+rhs.y + lhs.z*rhs.z;
}

template <typename T>
Vector3<T> abs() {
    return Vector3(std::abs(x), std::abs(y), std::abs(z));
}

template <typename T>
Vector3<T> Cross(const Vector3<T> &lhs, const Vector3<T> &rhs) {
    //convert to floats before subtraction to prevent floating point rounding errors
    double lhsx = lhs.x, lhsy = lhs.y, lhsz = lhs.z, rhsx = rhs.x, rhsy = rhs.y, rhsz = rhs.z;
    return Vector3(lhsy*rhsz - lhsz*rhsy, -(lhsx*rhsz - lhsz*rhsx), lhsx*rhsy - lhsy*rhsx);
}

template <typename T>
float Vector3<T>::Mag() const {
    return std::sqrt(x*x + y*y + z*z);
}

template <typename T>
Vector3<T> Unit(const Vector3<T> &vec) {
    return vec/vec.Mag();
}