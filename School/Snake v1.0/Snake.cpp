#include "Snake.h"


Snake::Snake(int x,int y)  {
 tail.resize(0);
 tail.push_back(TPoint(x,y));
 dir = right;
 food.resize(0);
}

void Snake::setDerection(Derection dir)  {
 this->dir = dir;
}

bool Snake::step() {
  bool add = false;
  if(food.size()) {
  add = food[food.size()-1].x == tail[tail.size()-1].x &&
  food[food.size()-1].y == tail[tail.size()-1].y;
  }
 for(int i = tail.size() - 1; i > 0;i--) {
  tail[i].x = tail[i-1].x;
  tail[i].y = tail[i-1].y;
 }
 if(dir == left || dir == right)
   tail[0].x += (dir == right)?1:-1;
 if(dir == top || dir == bottom)
   tail[0].y += (dir == bottom)?1:-1;
   if(add) {
    tail.push_back(food[0]);
    for(int i = 0;i< food.size();i++)  {
    food[i] = food[i +1];
    }
    food.pop_back();
   }
   for(int i = 1;i < tail.size();i++)  {
    if(tail[0].x == tail[i].x && tail[0].y == tail[i].y)
      return false;
   }
}


vector<TPoint> Snake::gettail()  {
return tail;
}

TPoint Snake::gethead()  {
 return tail[0];
}

Derection Snake::getDerection()  {
 return dir;
}

void Snake::eat(TPoint p) {
  food.push_back(p);

}

void Snake::reset(int x,int y) {
 tail.resize(0);
 tail.push_back(TPoint(x,y));
 dir = right;
 food.resize(0);
}
