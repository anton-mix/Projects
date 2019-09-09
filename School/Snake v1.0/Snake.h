#ifndef SNAKE_H
#define SNAKE_H
#include <vector>
#include <vcl.h>
using namespace std;

enum Derection {top,left,right,bottom};


class Snake  {
  Derection dir;
  vector<TPoint> tail;
  vector<TPoint> food;

  public:
  Snake(int x,int y);
  Snake();
  bool step();
  void setDerection(Derection);
  void eat(TPoint);
  void reset(int,int);
  Derection getDerection();
  vector<TPoint> gettail();
  TPoint gethead();
};

#endif