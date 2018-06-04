#ifndef MAGNET_H
#define MAGNET_H

#include <cmath>

class Stat {
public:
  Stat(){};
  void push_m(const double);
  double get_magnet() const;
  double get_susc() const;
  double get_binder() const;
  void reset();

private:
  unsigned t = 0;
  double m = 0;
  double m2 = 0;
  double m4 = 0;
};

void Stat::push_m(const double val) {
  m += fabs(val);
  double v2 = val * val;
  m2 += v2;
  m4 += v2 * v2;
  ++t;
}

double Stat::get_magnet() const { return m / t; }
double Stat::get_susc() const { return m2 / t - (m * m) / (t * t); }
double Stat::get_binder() const { return -m4 / t / (3 * m2 / t * m2 / t) + 1; }

void Stat::reset() {
  m = 0;
  m2 = 0;
  m4 = 0;
  t = 0;
}

#endif // MAGNET_H
