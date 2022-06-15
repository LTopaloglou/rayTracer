//#include "geometry.h"
//This is an implementation file because template classes cannot be implemented 
//in the .cpp and i want to keep the .h file clean

//Vector3 implementation
template <typename T> 
Vector3<T>::Vector3() {
    x = y = z = 0;
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