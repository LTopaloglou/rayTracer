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
    if (i == 0) return x;
    if (i == 1) return y;
    return z;
}

template <typename T>
T& Vector3<T>::operator[] (int i) {
    //check for invalid i, throw error if true
    Assert(i >= 0 && i <= 2);
    if (i == 0) return x;
    if (i == 1) return y;
    return z;
}

template <typename T>
Vector3<T> Vector3<T>::operator+ (const Vector3<T> &rhs) const {
    return Vector3<T>(x + rhs.x, y + rhs.y, z + rhs.z);
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
    return Vector3<T>(-x, -y, -z);
}

template <typename T>
Vector3<T> Vector3<T>::operator-(const Vector3<T> &rhs) const {
    return Vector3<T>(x-rhs.x, y-rhs.y, z-rhs.z);
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
    return Vector3<T>(x*i, y*i, z*i);
}

template <typename T>
Vector3<T>& Vector3<T>::operator*= (T i) {
    x *= i;
    y *= i;
    z *= i;
    return *this;
}

template <typename T> 
Vector3<T> operator* (T i, const Vector3<T> &vec) {
    return Vector3<T>(vec.x*i, vec.y*i, vec.z*i);
}

template <typename T>
Vector3<T>  Vector3<T>::operator/ (T i)const {
    //make sure were not dividing by zero
    Assert(i != 0);
    //calculate reciprocal of i first since this means only 1 division compared to 3
    //because division takes longer than multiplication
    float reciprocal = 1/i;
    return Vector3<T>(x*reciprocal, y*reciprocal, z*reciprocal);
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
Vector3<T> Vector3<T>::abs() const{
    return Vector3<T>(std::abs(x), std::abs(y), std::abs(z));
}

template <typename T>
Vector3<T> Cross(const Vector3<T> &lhs, const Vector3<T> &rhs) {
    //convert to floats before subtraction to prevent floating point rounding errors
    double lhsx = lhs.x, lhsy = lhs.y, lhsz = lhs.z, rhsx = rhs.x, rhsy = rhs.y, rhsz = rhs.z;
    return Vector3<T>(lhsy*rhsz - lhsz*rhsy, -(lhsx*rhsz - lhsz*rhsx), lhsx*rhsy - lhsy*rhsx);
}

template <typename T>
float Vector3<T>::Mag() const {
    return std::sqrt(x*x + y*y + z*z);
}

template <typename T>
Vector3<T> Unit(const Vector3<T> &vec) {
    return vec/vec.Mag();
}

template <typename T>
Point3<T>::Point3() {
    x = y = z = 0;
}

template <typename T>
Point3<T>::Point3(T x_, T y_, T z_) {
    x = x_;
    y = y_;
    z = z_;
}

template <typename T>
T Point3<T>::operator[](int i) const {
    //check for invalid i, throw error if true
    Assert(i >= 0 && i <= 2);
    if (i == 0) return x;
    if (i == 1) return y;
    return z;
}

template <typename T>
T& Point3<T>::operator[] (int i) {
    //check for invalid i, throw error if true
    Assert(i >= 0 && i <= 2);
    if (i == 0) return x;
    if (i == 1) return y;
    return z;
}

template <typename T>
Point3<T> Point3<T>::operator+(const Vector3<T> &vec) const {
    return Point3(x + vec.x, y + vec.y, z + vec.z);
} 

template <typename T>
Point3<T>& Point3<T>::operator+=(const Vector3<T> &vec) {
    x += vec.x;
    y += vec.y;
    z += vec.z;
    return *this;
}

template <typename T>
Point3<T> Point3<T>::operator-(const Vector3<T> &vec) const {
    return Point3(x - vec.x, y - vec.y, z = vec.z);
}

template <typename T>
Point3<T>& Point3<T>::operator-=(const Vector3<T> &vec) {
    x -= vec.x;
    y -= vec.y;
    z -= vec.z;
    return *this;
}

template <typename T>
Vector3<T> Point3<T>::operator-(const Point3<T> &rhs) const {
    return Vector3<T>(x - rhs.x, y - rhs.y, z - rhs.z);
}

template <typename T>
float Distance(const Point3<T> &lhs, const Point3<T> &rhs) {
    return (lhs-rhs).Mag();
}

template <typename T>
Point3<float> Lerp(const Point3<T> &lhs, const Point3<T> &rhs, float n) {
    //n values between 0 and 1 interpolate, negative and >1 values extrapolate
    Vector3<float> vec = rhs - lhs;
    return lhs + (n*vec);
}

template <typename T>
Normal3<T>::Normal3(const Vector3<T> & vec) {
    x = vec.x;
    y = vec.y;
    z = vec.z;
}

Ray::Ray() {
    max = std::numeric_limits<float>::max();
    time = 0.f;
    //medium = nullptr; re implement when mediums work
}

Ray::Ray(const Point3f & origin_, const Vector3f & direction_) {
    origin = origin_;
    direction = direction_; 
    max = std::numeric_limits<float>::max();
    time = 0.f;
    //medium = nullptr; re implement when mediums work
}

Point3f Ray::operator()(float f) {
    return origin + (direction * f);
}

template <typename T>
Bounds3<T>::Bounds3() {
    //default constructor sets points to invalid
    //values since max < min
    T lowest = std::numeric_limits<T>::min();
    T highest = std::numeric_limits<T>::max();
    min = Point3<T>(highest, highest, highest);
    max = Point3<T>(lowest, lowest, lowest);
}

template <typename T>
Bounds3<T>::Bounds3(const Point3<T> & min_, T x_, T y_, T z_) {
    min = min_;
    max = Point3<T>(min.x + x_, min.y + y_, min.z + z_);
}

template <typename T>
Bounds3<T>::Bounds3(const Point3<T> & p1, const Point3<T> & p2) {
    min = Point3<T> (std::min(p1.x, p2.x), std::min(p1.y, p2.y), std::min(p1.z, p2.z));
    max = Point3<T> (std::max(p1.x, p2.x), std::max(p1.y, p2.y), std::max(p1.z, p2.z));
}

template <typename T>
Bounds3<T>::Bounds3(const Point3<T> & point) {
    min = max = point;
}

Matrix4x4::Matrix4x4() {
    //create identity matrix
    m[0][0] = m[1][1] = m[2][2] = m[3][3] = 1.f;
    m[0][1] = m[0][2] = m[0][3] = 0.f;
    m[1][0] = m[1][2] = m[1][3] = 0.f;
    m[2][0] = m[2][1] = m[2][3] = 0.f;
    m[3][0] = m[3][1] = m[3][2] = 0.f;
}

Matrix4x4::Matrix4x4(float m00, float m01, float m02, float m03,
                    float m10, float m11, float m12, float m13,
                    float m20, float m21, float m22, float m23,
                    float m30, float m31, float m32, float m33) {
    m[0][0] = m00;
    m[0][1] = m01;
    m[0][2] = m02;
    m[0][3] = m03;

    m[1][0] = m10;
    m[1][1] = m11;
    m[1][2] = m12;
    m[1][3] = m13;

    m[2][0] = m20;
    m[2][1] = m21;
    m[2][2] = m22;
    m[2][3] = m23;

    m[3][0] = m30;
    m[3][1] = m31;
    m[3][2] = m32;
    m[3][3] = m33;
}

Matrix4x4 Matrix4x4::operator+(const Matrix4x4 & rhs) const {
    return Matrix4x4(m[0][0] + rhs.m[0][0], m[0][1] + rhs.m[0][1], m[0][2] + rhs.m[0][2], m[0][3] + rhs.m[0][3],
                     m[1][0] + rhs.m[1][0], m[1][1] + rhs.m[1][1], m[1][2] + rhs.m[1][2], m[1][3] + rhs.m[1][3],
                     m[2][0] + rhs.m[2][0], m[2][1] + rhs.m[2][1], m[2][2] + rhs.m[2][2], m[2][3] + rhs.m[2][3],
                     m[3][0] + rhs.m[3][0], m[3][1] + rhs.m[3][1], m[3][2] + rhs.m[3][2], m[3][3] + rhs.m[3][3]);
}

Matrix4x4 Matrix4x4::operator+=(const Matrix4x4 & rhs) {
    m[0][0] += rhs.m[0][0]; m[0][1] += rhs.m[0][1]; m[0][2] += rhs.m[0][2]; m[0][3] += rhs.m[0][3];
    m[1][0] += rhs.m[1][0]; m[1][1] += rhs.m[1][1]; m[1][2] += rhs.m[1][2]; m[1][3] += rhs.m[1][3];
    m[2][0] += rhs.m[2][0]; m[2][1] += rhs.m[2][1]; m[2][2] += rhs.m[2][2]; m[2][3] += rhs.m[2][3];
    m[3][0] += rhs.m[3][0]; m[3][1] += rhs.m[3][1]; m[3][2] += rhs.m[3][2]; m[3][3] += rhs.m[3][3];
    return *this; 
}

Matrix4x4 Matrix4x4::operator-(const Matrix4x4 & rhs) const {
    return Matrix4x4(m[0][0] - rhs.m[0][0], m[0][1] - rhs.m[0][1], m[0][2] - rhs.m[0][2], m[0][3] - rhs.m[0][3],
                     m[1][0] - rhs.m[1][0], m[1][1] - rhs.m[1][1], m[1][2] - rhs.m[1][2], m[1][3] - rhs.m[1][3],
                     m[2][0] - rhs.m[2][0], m[2][1] - rhs.m[2][1], m[2][2] - rhs.m[2][2], m[2][3] - rhs.m[2][3],
                     m[3][0] - rhs.m[3][0], m[3][1] - rhs.m[3][1], m[3][2] - rhs.m[3][2], m[3][3] - rhs.m[3][3]);
}

Matrix4x4 Matrix4x4::operator-=(const Matrix4x4 & rhs) {
    m[0][0] -= rhs.m[0][0]; m[0][1] -= rhs.m[0][1]; m[0][2] -= rhs.m[0][2]; m[0][3] -= rhs.m[0][3];
    m[1][0] -= rhs.m[1][0]; m[1][1] -= rhs.m[1][1]; m[1][2] -= rhs.m[1][2]; m[1][3] -= rhs.m[1][3];
    m[2][0] -= rhs.m[2][0]; m[2][1] -= rhs.m[2][1]; m[2][2] -= rhs.m[2][2]; m[2][3] -= rhs.m[2][3];
    m[3][0] -= rhs.m[3][0]; m[3][1] -= rhs.m[3][1]; m[3][2] -= rhs.m[3][2]; m[3][3] -= rhs.m[3][3];
    return *this; 
}

template <typename T>
Vector3f Matrix4x4::operator*(const Vector3<T> & rhs) const {
    //treat matrix like 3x3 matrix for this operation
    return Vector3f(m[0][0]*rhs.x + m[0][1]*rhs.y + m[0][2]*rhs.z,
                    m[1][0]*rhs.x + m[1][1]*rhs.y + m[1][2]*rhs.z,
                    m[2][0]*rhs.x + m[2][1]*rhs.y + m[2][2]*rhs.z);
}

Matrix4x4 Matrix4x4::operator*(const Matrix4x4 & rhs) const {
    //multiply rows of this with columns of rhs for each element
    return Matrix4x4(m[0][0]*rhs.m[0][0] + m[0][1]*rhs.m[1][0] + m[0][2]*rhs.m[2][0] + m[0][3]*rhs.m[3][0], m[0][0]*rhs.m[0][1] + m[0][1]*rhs.m[1][1] + m[0][2]*rhs.m[2][1] + m[0][3]*rhs.m[3][1], m[0][0]*rhs.m[0][2] + m[0][1]*rhs.m[1][2] + m[0][2]*rhs.m[2][2] + m[0][3]*rhs.m[3][2], m[0][0]*rhs.m[0][3] + m[0][1]*rhs.m[1][3] + m[0][2]*rhs.m[2][3] + m[0][3]*rhs.m[3][3],
                     m[1][0]*rhs.m[0][0] + m[1][1]*rhs.m[1][0] + m[1][2]*rhs.m[2][0] + m[1][3]*rhs.m[3][0], m[1][0]*rhs.m[0][1] + m[1][1]*rhs.m[1][1] + m[1][2]*rhs.m[2][1] + m[1][3]*rhs.m[3][1], m[1][0]*rhs.m[0][2] + m[1][1]*rhs.m[1][2] + m[1][2]*rhs.m[2][2] + m[1][3]*rhs.m[3][2], m[1][0]*rhs.m[0][3] + m[1][1]*rhs.m[1][3] + m[1][2]*rhs.m[2][3] + m[1][3]*rhs.m[3][3],
                     m[2][0]*rhs.m[0][0] + m[2][1]*rhs.m[1][0] + m[2][2]*rhs.m[2][0] + m[2][3]*rhs.m[3][0], m[2][0]*rhs.m[0][1] + m[2][1]*rhs.m[1][1] + m[2][2]*rhs.m[2][1] + m[2][3]*rhs.m[3][1], m[2][0]*rhs.m[0][2] + m[2][1]*rhs.m[1][2] + m[2][2]*rhs.m[2][2] + m[2][3]*rhs.m[3][2], m[2][0]*rhs.m[0][3] + m[2][1]*rhs.m[1][3] + m[2][2]*rhs.m[2][3] + m[2][3]*rhs.m[3][3],
                     m[3][0]*rhs.m[0][0] + m[3][1]*rhs.m[1][0] + m[3][2]*rhs.m[2][0] + m[3][3]*rhs.m[3][0], m[3][0]*rhs.m[0][1] + m[3][1]*rhs.m[1][1] + m[3][2]*rhs.m[2][1] + m[3][3]*rhs.m[3][1], m[3][0]*rhs.m[0][2] + m[3][1]*rhs.m[1][2] + m[3][2]*rhs.m[2][2] + m[3][3]*rhs.m[3][2], m[3][0]*rhs.m[0][3] + m[3][1]*rhs.m[1][3] + m[3][2]*rhs.m[2][3] + m[3][3]*rhs.m[3][3]);
}

Matrix4x4 Matrix4x4::operator*(float s) const {
    return Matrix4x4(m[0][0]*s, m[0][1]*s, m[0][2]*s, m[0][3]*s,
                     m[1][0]*s, m[1][1]*s, m[1][2]*s, m[1][3]*s,
                     m[2][0]*s, m[2][1]*s, m[2][2]*s, m[2][3]*s,
                     m[3][0]*s, m[3][1]*s, m[3][2]*s, m[3][3]*s);
}

Matrix4x4 Matrix4x4::operator*= (float s) {
    m[0][0]*=s; m[0][1]*=s; m[0][2]*=s; m[0][3]*=s;
    m[1][0]*=s; m[1][1]*=s; m[1][2]*=s; m[1][3]*=s;
    m[2][0]*=s; m[2][1]*=s; m[2][2]*=s; m[2][3]*=s;
    m[3][0]*=s; m[3][1]*=s; m[3][2]*=s; m[3][3]*=s;
    return *this;
}

void Matrix4x4::print() const {
    std::cout << m[0][0] << " " << m[0][1] << " " << m[0][2] << " " << m[0][3] << std::endl;
    std::cout << m[1][0] << " " << m[1][1] << " " << m[1][2] << " " << m[1][3] << std::endl;
    std::cout << m[2][0] << " " << m[2][1] << " " << m[2][2] << " " << m[2][3] << std::endl;
    std::cout << m[3][0] << " " << m[3][1] << " " << m[3][2] << " " << m[3][3] << std::endl;
}

Matrix4x4 Matrix4x4::transpose() const {
    return Matrix4x4(m[0][0], m[1][0], m[2][0], m[3][0],
                     m[0][1], m[1][1], m[2][1], m[3][1],
                     m[0][2], m[1][2], m[2][2], m[3][2],
                     m[0][3], m[1][3], m[2][3], m[3][3]);
}

Matrix4x4 Matrix4x4::inverse() const {
    return this->adjugate()*this->determinant();
}

Matrix4x4 Matrix4x4::adjugate() const {
    return Matrix4x4(m[1][1]*m[2][2]*m[3][3] + m[2][1]*m[3][2]*m[1][3] + m[3][1]*m[1][2]*m[2][3] - m[3][1]*m[2][2]*m[1][3] - m[2][1]*m[1][2]*m[3][3] - m[1][1]*m[3][2]*m[2][3],
                   -(m[1][0]*m[2][2]*m[3][3] + m[2][0]*m[3][2]*m[1][3] + m[3][0]*m[1][2]*m[2][3] - m[3][0]*m[2][2]*m[1][3] - m[2][0]*m[1][2]*m[3][3] - m[1][0]*m[3][2]*m[2][3]),
                     m[1][0]*m[2][1]*m[3][3] + m[2][0]*m[3][1]*m[1][3] + m[3][0]*m[1][1]*m[2][3] - m[3][0]*m[2][1]*m[1][3] - m[2][0]*m[1][1]*m[3][3] - m[1][0]*m[3][1]*m[2][3],
                   -(m[1][0]*m[2][1]*m[3][2] + m[2][0]*m[3][1]*m[1][2] + m[3][0]*m[1][1]*m[2][2] - m[3][0]*m[2][1]*m[1][2] - m[2][0]*m[1][1]*m[3][2] - m[1][0]*m[3][1]*m[2][2]),
                     
                   -(m[0][1]*m[2][2]*m[3][3] + m[2][1]*m[3][2]*m[0][3] + m[3][1]*m[0][2]*m[2][3] - m[3][1]*m[2][2]*m[0][3] - m[2][1]*m[0][2]*m[3][3] - m[0][1]*m[3][2]*m[2][3]),
                     m[0][0]*m[2][2]*m[3][3] + m[2][0]*m[3][2]*m[0][3] + m[3][0]*m[0][2]*m[2][3] - m[3][0]*m[2][2]*m[0][3] - m[2][0]*m[0][2]*m[3][3] - m[0][0]*m[3][2]*m[2][3],
                   -(m[0][0]*m[2][1]*m[3][3] + m[2][0]*m[3][1]*m[0][3] + m[3][0]*m[0][1]*m[2][3] - m[3][0]*m[2][1]*m[0][3] - m[2][0]*m[0][1]*m[3][3] - m[0][0]*m[3][1]*m[2][3]),
                     m[0][0]*m[2][1]*m[3][2] + m[2][0]*m[3][1]*m[0][2] + m[3][0]*m[0][1]*m[2][2] - m[3][0]*m[2][1]*m[0][2] - m[2][0]*m[0][1]*m[3][2] - m[0][0]*m[3][1]*m[2][2],

                     m[0][1]*m[1][2]*m[3][3] + m[1][1]*m[3][2]*m[0][3] + m[3][1]*m[0][2]*m[1][3] - m[3][1]*m[1][2]*m[0][3] - m[1][1]*m[0][2]*m[3][3] - m[0][1]*m[3][2]*m[1][3],
                   -(m[0][0]*m[1][2]*m[3][3] + m[1][0]*m[3][2]*m[0][3] + m[3][0]*m[0][2]*m[1][3] - m[3][0]*m[1][2]*m[0][3] - m[1][0]*m[0][2]*m[3][3] - m[0][0]*m[3][2]*m[1][3]),
                     m[0][0]*m[1][1]*m[3][3] + m[1][0]*m[3][1]*m[0][3] + m[3][0]*m[0][1]*m[1][3] - m[3][0]*m[1][1]*m[0][3] - m[1][0]*m[0][1]*m[3][3] - m[0][0]*m[3][1]*m[1][3],
                   -(m[0][0]*m[1][1]*m[3][2] + m[1][0]*m[3][1]*m[0][2] + m[3][0]*m[0][1]*m[1][2] - m[3][0]*m[1][1]*m[0][2] - m[1][0]*m[0][1]*m[3][2] - m[0][0]*m[3][1]*m[1][2]),

                   -(m[0][1]*m[1][2]*m[2][3] + m[1][1]*m[2][2]*m[0][3] + m[2][1]*m[0][2]*m[1][3] - m[2][1]*m[1][2]*m[0][3] - m[1][1]*m[0][2]*m[2][3] - m[0][1]*m[2][2]*m[1][3]),
                     m[0][0]*m[1][2]*m[2][3] + m[1][0]*m[2][2]*m[0][3] + m[2][0]*m[0][2]*m[1][3] - m[2][0]*m[1][2]*m[0][3] - m[1][0]*m[0][2]*m[2][3] - m[0][0]*m[2][2]*m[1][3],
                   -(m[0][0]*m[1][1]*m[2][3] + m[1][0]*m[2][1]*m[0][3] + m[2][0]*m[0][1]*m[1][3] - m[2][0]*m[1][1]*m[0][3] - m[1][0]*m[0][1]*m[2][3] - m[0][0]*m[2][1]*m[1][3]),
                     m[0][0]*m[1][1]*m[2][2] + m[1][0]*m[2][1]*m[0][2] + m[2][0]*m[0][1]*m[1][2] - m[2][0]*m[1][1]*m[0][2] - m[1][0]*m[0][1]*m[2][2] - m[0][0]*m[2][1]*m[1][2]
                    );
}

float Matrix4x4::determinant() const {
    return ((m[0][0]*m[1][1]*m[2][2]*m[3][3]) + (m[0][0]*m[2][1]*m[3][2]*m[1][3]) + (m[0][0]*m[3][1]*m[1][2]*m[2][3])
           -(m[0][0]*m[3][1]*m[2][2]*m[1][3]) - (m[0][0]*m[2][1]*m[1][2]*m[3][3]) - (m[0][0]*m[1][1]*m[3][2]*m[2][3])
           -(m[1][0]*m[0][1]*m[2][2]*m[3][3]) - (m[2][0]*m[0][1]*m[3][2]*m[1][3]) - (m[3][0]*m[0][1]*m[1][2]*m[2][3])
           +(m[3][0]*m[0][1]*m[2][2]*m[1][3]) + (m[2][0]*m[0][1]*m[1][2]*m[3][3]) + (m[1][0]*m[0][1]*m[3][2]*m[2][3])
           +(m[1][0]*m[2][1]*m[0][2]*m[3][3]) + (m[2][0]*m[3][1]*m[0][2]*m[1][3]) + (m[3][0]*m[1][1]*m[0][2]*m[2][3])
           -(m[3][0]*m[2][1]*m[0][2]*m[1][3]) - (m[2][0]*m[1][1]*m[0][2]*m[3][3]) - (m[1][0]*m[3][1]*m[0][2]*m[2][3])
           -(m[1][0]*m[2][1]*m[3][2]*m[0][3]) - (m[2][0]*m[3][1]*m[1][2]*m[0][3]) - (m[3][0]*m[1][1]*m[2][2]*m[0][3])
           +(m[3][0]*m[2][1]*m[1][2]*m[0][3]) + (m[2][0]*m[1][1]*m[3][2]*m[0][3]) + (m[1][0]*m[3][1]*m[2][2]*m[0][3]));
}

template <typename T>
Matrix4x4 operator* (T s, const Matrix4x4 &mat) {
    return mat*s;
}