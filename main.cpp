//============================================================
// file: main.cpp
//============================================================
#include <iostream>
#include <cstring>
#include <initializer_list>
#include <cassert>

//MATRIX and VECTOR classes assignment
#include "vector3d_T.h"
#include "matrix3d_T.h"


template <typename T>
void print(T v) { std::cout << v << "\n"; }

template <typename T>
void show_vect(T v) { std::cout << v.name() << " is: " << v << "\n"; }

template <typename T>
void show_mat(T m) { std::cout << m.name() << " is: " << m << "\n"; }

void test_vectors() {
  print("\n====================  TESTING VECTORS  ========================");
  vector3d<double> u("u", 3, {1,  2,  4});
  std::cout << u.name() << "\n";
  std::cout << u << "\n";
  u.zero();
  show_vect(u);
  vector3D v("v", 3, {8, 16, 32});
  vector3D i("i", 3, {1, 0, 0}), j("j", 3, {0, 1, 0}), k("k", 3, {0, 0, 1});
  vector3D w(3 * i + 4 * j - 2 * k);

  show_vect(u);
  show_vect(v);
  show_vect(i);
  show_vect(j);
  show_vect(k);
  j + k;
  show_vect(w);

  assert(u == u);
  assert(u != v);
  
  assert(u + v == v + u);
  assert(u - v == -(v - u));
  assert(-(-u) == u);
  
  assert(3.0 + u == u + 3.0);
  assert(3.0 * u == u * 3.0);
  assert((u - 3.0) == -(3.0 - u));
  assert((5.0 * u) / 5.0 == u);

  assert(u + vector3D::zero() == u);

  assert(i.dot(j) == j.dot(k) == k.dot(i) == 0);
  
  assert(i.cross(j) == k);
  assert(j.cross(k) == i);
  assert(k.cross(i) == j);
  
  assert(u.cross(v) == -v.cross(u));
  assert(u.cross(v + w) == u.cross(v) + u.cross(w));
  assert((u.cross(v)).dot(u) == 0);

  print(i.angle(j));
  print(M_PI/2);

  assert(i.angle(j) == M_PI_2);
  assert(j.angle(k) == M_PI_2);
  assert(k.angle(i) == M_PI_2);
  
  vector3D uhat = u / u.magnitude();
  show_vect(u);
  show_vect(uhat);
  print(uhat.magnitude());
  assert(uhat.magnitude() - 1.0 < 1.0e-10);

  print("...test vectors assertions passed");
  print("====================  FINISHED testing vectors  ========================");
}

void test_matrices() {
  print("\n====================  TESTING MATRICES  ========================");
  matrix3dD a("a", 3, {3, 2, 0,   0, 0, 1,   2, -2, 1});
  matrix3dD b("b", 3, {1, 0, 5,   2, 1, 6,   3,  4, 0});
  
  matrix3dD id = matrix3dD::identity(3);
  assert(a * id == a);
  assert(a * b != -b * a);
  assert((a * b).transpose() == b.transpose() * a.transpose());

  matrix3dD acopy(a);    // copy constructor
  matrix3dD a2copy = a;  // copy constructor

  matrix3dD bcopy;
  bcopy = b;        // assignment operator

  matrix3dD ainv = a.inverse();
  matrix3dD binv = b.inverse();

  show_mat(a);
  show_mat(b);
  show_mat(-a);
  show_mat(-b);
  show_mat(a * b);
  printf("|a| = %.2f\n", a.determinant());
  printf("|b| = %.2f\n", b.determinant());
  show_mat(a.transpose());
  show_mat(b.transpose());

  show_mat(a.minors());
  show_mat(b.minors());

  show_mat(a.cofactor());
  show_mat(b.cofactor());

  show_mat(a.adjugate());
  show_mat(b.adjugate());


  show_mat(ainv);
  show_mat(binv);
  show_mat(a * ainv);
  show_mat(b * binv);
  show_mat(matrix3dD::identity(3));

  assert(a * ainv == matrix3dD::identity(3));
  assert(a * ainv == ainv * a);
  assert(b * binv == matrix3dD::identity(3));
  assert(b * binv == binv * b);
  assert(a.transpose().transpose() == a);
  assert(a.determinant() == a.transpose().determinant());

  assert(a + b == b + a);
  assert(a - b == -(b - a));
  assert(3.0 + a == a + 3.0);
  assert(3.0 * a == a * 3.0);
  assert((a + 3.0) - 3.0 == a);
  assert((3.0 * a) / 3.0 == a);
  assert(-(-a) == a);

  matrix3dD zerod("zerod", 3, {1, 2, 3,   4, 5, 6,   7, 8, 9});
  assert(zerod.determinant() == 0);

  print("...test matrices assertions passed");
  print("====================  FINISHED testing matrices  ========================");
}

void test_matrices_and_vectors() {
  print("\n====================  TESTING MATRICES and VECTORS  ========================");
  vector3D p("p", 2, {1, 2});
  matrix3dD m("m", 2, {1, 2,   3, 4});
  show_vect(p);
  show_mat(m);
  assert(p * m == m * p);

  vector3D q("q", 3, {1, 2, 3});
  matrix3dD n("n", 3, {1, 2, 3,   4, 5, 6,   7, 8, 9});
  show_vect(q);
  show_mat(n);
  assert(q * n == n * q);
  print("...test_matrices_and_vectors assertions passed");
  print("====================  FINISHED testing matrices and vectors  ========================");
}


int main(int argc, const char * argv[]) {
  test_vectors();
  test_matrices();
  test_matrices_and_vectors();
    
  return 0;
}


