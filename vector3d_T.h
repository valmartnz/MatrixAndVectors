//============================================================
// FILE: vector_3dT.h
//============================================================

#ifndef __vector3d_T_H__
#define __vector3d_T_H__

#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <initializer_list>

template <typename T> class vector3d;
template <typename T> std::ostream& operator<<(std::ostream& os, const vector3d<T>& v);

typedef vector3d<double> vector3D;
typedef vector3d<float>  vector3F;
typedef vector3d<int>    vector3I;
typedef vector3d<long>   vector3L;
typedef vector3D vec3;


template <typename T>
class vector3d {  // class that serves as both vectors and points
public:
  vector3d();
  vector3d(const std::string& name, int dims);
  vector3d(const std::string& name, int dims, const std::initializer_list<T>& li);
//---------------------------------------------------------------------
  T  operator[](int i) const;
  T& operator[](int i);
//---------------------------------------------------------------------
  std::string name() const;
  void name(const std::string& name);

  double x() const { return data_[0]; }
  double y() const  { return data_[1]; }
  double z() const  { return data_[2]; }

  double& x() { return data_[0]; }
  double& y() { return data_[1]; }
  double& z() { return data_[2]; }

//---------------------------------------
  vector3d<T>& operator+=(const vector3d<T>& v);
  vector3d<T>& operator-=(const vector3d<T>& v);
//---------------------------------------------------------------------
  vector3d<T>& operator+=(T k);
  vector3d<T>& operator-=(T k);
  vector3d<T>& operator*=(T k);
  vector3d<T>& operator/=(T k);
//---------------------------------------------------------------------
  vector3d<T> operator-();
//---------------------------------------------------------------------
  friend vector3d<T> operator+(const vector3d<T>& u, const vector3d<T>& v) {
    check_equal_dims(u, v);
    return vector3d<T>(u.name_ + "+" + v.name_, u.dims_,
                    { u[0] + v[0], u[1] + v[1], u[2] + v[2], 0} );
  }  

  friend vector3d<T> operator-(const vector3d<T>& u, const vector3d<T>& v) { /* TODO */ }
//---------------------------------------------------------------------
  friend vector3d<T> operator+(T k, const vector3d<T>& v) {
    return vector3d<T>(std::to_string(k) + "+" + v.name_, v.dims_, { k + v[0], k + v[1], k + v[2], 0 });
  }
  friend vector3d<T> operator+(const vector3d<T>& v, T k) { /* TODO */ }
//---------------------------------------------------------------------
  friend vector3d<T> operator-(T k, const vector3d<T>& v) { /* TODO */ }
  friend vector3d<T> operator-(const vector3d<T>& v, T k) { /* TODO */ }
  //---------------------------------------------------------------------
  friend vector3d<T> operator*(T k, const vector3d<T>& v) {
    return vector3d<T>(std::to_string(k) + v.name_, v.dims_, { k * v[0], k * v[1], k * v[2], 0 });
  }
  friend vector3d<T> operator*(const vector3d<T>& v, T k) { /* TODO */ };
  //---------------------------------------------------------------------
  friend vector3d<T> operator/(const vector3d<T>& v, T k) {
    if (k == 0) { throw new std::invalid_argument("divide by zero"); }
    double kinv = 1.0 / k;
    return kinv * v;
  }
//---------------------------------------------------------------------
  friend bool operator<(const vector3d<T>& u, const vector3d<T>& v) {
    check_equal_dims(u, v);
    return u.magnitude() < v.magnitude();
  }
  friend bool operator==(const vector3d<T>& u, const vector3d<T>& v) {
    check_equal_dims(u, v);
    return /* u.name_ == v.name_ && */ u[0] == v[0] && u[1] == v[1] && u[2] == v[2];
  }
  friend bool operator!=(const vector3d<T>& u, const vector3d<T>& v) { return !(u == v); }
//---------------------------------------------------------------------
  double dot(const vector3d<T>& other) const;
  double magnitude() const;
  double norm() const { return magnitude(); }  // L2 norm
  double angle(const vector3d<T>& other) const;
  vector3d<T> cross(const vector3d<T>& other) const;
//---------------------------------------------------------------------
  static vector3d<T> zero();
//---------------------------------------------------------------------
  friend std::ostream& operator<< <>(std::ostream& os, const vector3d<T>& v);

private:
  friend void check_equal_dims(const vector3d& u, const vector3d& v) {
    if (u.dims_ != v.dims_) { throw new std::invalid_argument("vector3d dims mismatch"); }
  }
  void check_bounds(int i) const;

private:
  std::string name_;
  int dims_;
  double data_[4];
};


//==============================================================================

template <typename T> vector3d<T>::vector3d() : vector3d("no_name", 3) {}  // 3d default dims
template <typename T> vector3d<T>::vector3d(const std::string& name, int dims)
  : name_(name), dims_(dims) {
    memset(data_, 0, dims_ * sizeof(double));
    data_[3] = 0;  // vectors have 0 at end, pts have 1
  }
template <typename T> vector3d<T>::vector3d(const std::string& name, int dims,
                                            const std::initializer_list<T>& li)
  : vector3d(name, dims) {
    int i = 0;
    for (T value : li) {
      if (i > dims_) { break; }
      data_[i++] = value;
    }
    data_[3] = 0;
  }

template <typename T> std::string vector3d<T>::name() const { return name_; }
template <typename T> void vector3d<T>::name(const std::string& name) { name_ = name; }

//---------------------------------------------------------------------
double x() const { return data_[0]; }
double y() const { /* TODO */ }
double z() const { /* TODO */ }

double& x() { return data_[0]; }
double& y() { /* TODO */ }
double& z() { /* TODO */ }

//---------------------------------------------------------------------

template <typename T>
vector3d<T>& vector3d<T>::operator+=(const vector3d<T>& v) {
  vector3d<T>& u = *this;
  u[0] += v[0];  u[1] += v[1];  u[2] += v[2];  u[3] = v[3];
  return *this;
}
template <typename T>
vector3d<T>& vector3d<T>::operator-=(const vector3d<T>& v) { /* TODO */ }
//---------------------------------------------------------------------
template <typename T> vector3d<T>& vector3d<T>::operator+=(T k) { /* TODO */ }
template <typename T> vector3d<T>& vector3d<T>::operator-=(T k) { /* TODO */ }
template <typename T> vector3d<T>& vector3d<T>::operator*=(T k) { /* TODO */ }
template <typename T> vector3d<T>& vector3d<T>::operator/=(T k) { /* TODO */ }

//---------------------------------------------------------------------
template <typename T>  /* read only idx */
T  vector3d<T>::operator[](int i) const {  check_bounds(i);  return data_[i]; }

template <typename T> T& vector3d<T>::operator[](int i) { /* TODO */ } // rw idx

//-----------------------
template <typename T>
vector3d<T> vector3d<T>::operator-() { return vector3d("-" + name_, dims_, { -data_[0], -data_[1], -data_[2], 0 }); }
//-----------------------
template <typename T>
double vector3d<T>::dot(const vector3d<T>& v) const { /* TODO */ }
//-----------------------
template <typename T>
double vector3d<T>::magnitude() const {  return sqrt(dot(*this));  }

template <typename T>
double vector3d<T>::angle(const vector3d<T>& v) const { /* TODO */ }
//-----------------------
template <typename T>
vector3d<T> vector3d<T>::cross(const vector3d<T>& v) const {
  check_equal_dims(*this, v);
  if (v.dims_ != 3) {
    throw new std::invalid_argument("cross_product only implemented for vector3d's");
  }
  return vector3d<T>(name_ + " x " + v.name_, dims_,
                    { data_[1]*v[2] - data_[2]*v[1], -(data_[0]*v[2] - data_[2]*v[0]),
                      data_[0]*v[1] - data_[1]*v[0], 0 });
}
//-----------------------
template <typename T>
vector3d<T> vector3d<T>::zero() { return vector3d<T>("zero", 3, {0, 0, 0, 0}); }
//-----------------------

static double value(double val) { return abs(val) < 1e-5 ? 0 : val; }

template <typename T>
std::ostream& operator<<(std::ostream& os, const vector3d<T>& v) {
  os << "<" << v.name_ << ", ";
  if (v.dims_ == 0) { os << "empty>"; }
  else {
    for (int i = 0; i < v.dims_ + 1; ++i) {
      os << std::setw(3) << std::setprecision(3) << value(v[i]);
      if (i < v.dims_) { os << " "; }
    }
    os << ">";
  }
  return os;
}

template <typename T>
void vector3d<T>::check_bounds(int i) const {  // 1 extra dimension for pts/vectors
  if (i > dims_) { throw new std::invalid_argument("out of bounds"); }
}


#endif
