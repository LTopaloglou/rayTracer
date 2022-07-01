#include "geometry.h"

int main() {
    Matrix4x4 mat(0.f, 1.f, 2.f, 3.f,
                  4.f, 5.f, 6.f, 7.f,
                  8.f, 9.f, 10.f, 11.f,
                  12.f, 13.f, 14.f, 15.f);
    Matrix4x4 mat2 = mat.adjugate();
    //mat.print();
    Matrix4x4 identity;
    mat2.print();
    std::cout << "inverse:" << std::endl;
    Matrix4x4 mat3 = mat.inverse();
    mat3.print();
    return 0;
}