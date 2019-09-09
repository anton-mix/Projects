//---------------------------------------------------------------------------
#include <vcl.h>
#pragma hdrstop

#include <cstdlib>
#include <ctime>
#include <vector>
#include "xim.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TForm1 *Form1;
//---------------------------------------------------------------------------
__fastcall TForm1::TForm1(TComponent* Owner)
        : TForm(Owner)  {}
//---------------------------------------------------------------------------
int count = 0,dx = -1,dy = 0;
std::vector <TPoint> snake;
TPoint apple;
const int  width = 20, height = 20, x = 540,y = 340;
//---------------------------------------------------------------------------
void apprand() {
  srand(time(NULL));
  bool prov = true;
  while(prov)  {
    prov = false;
    apple.x = std::rand()%(x-40)/20;
    apple.y = std::rand()%(y-40)/20;
    for(int i = 0;i < snake.size();i ++)
      if(apple.x == snake[i].x && apple.y == snake[i].y)
        prov == true;
  }
}
//---------------------------------------------------------------------------
void __fastcall TForm1::FormCreate(TObject *Sender)  {
  snake.clear();
  snake.resize(3,TPoint(0,0));
  for(int i = 0;i<snake.size();i++)  {
    snake[i].x=11+i;
    snake[i].y=7;
    }
  apprand();
  Form1->DoubleBuffered = true;
}
//---------------------------------------------------------------------------
void __fastcall TForm1::FormResize(TObject *Sender)  {
  Form1->Invalidate();
}
//---------------------------------------------------------------------------
void __fastcall TForm1::FormPaint(TObject *Sender)  {
  Form1->Color = (TColor)RGB(222,222,222);

  Form1->Canvas->Brush->Style = bsSolid;
  Form1->Canvas->Pen->Style = psClear;
  Form1->Canvas->Brush->Color = (TColor)RGB(255,255,255);
  Form1->Canvas->Rectangle(20,20,x-20,y-20);

  Form1->Canvas->Pen->Style = psClear;
  Form1->Canvas->Font->Name = "Arial";
  Form1->Canvas->Font->Size = 80;
  Form1->Canvas->Font->Color = (TColor)RGB(222,222,222);
  Form1->Canvas->TextOutA(
  (Form1->ClientWidth - Form1->Canvas->TextWidth(count))/2,
  (Form1->ClientHeight - Form1->Canvas->TextHeight(count))/2,
  count);

  Form1->Canvas->Brush->Color = (TColor)RGB(255,106,106);
  Form1->Canvas->Ellipse((apple.x*20+30)-10,(apple.y*20+30)-10,(apple.x*20+30)+10,(apple.y*20+30)+10);

  Form1->Canvas->Brush->Color = (TColor)RGB(139,137,137);
  for(int i = 0;i < snake.size();i ++)
    Form1->Canvas->Ellipse((snake[i].x*20+30)-10,(snake[i].y*20+30)-10,(snake[i].x*20+30)+10,(snake[i].y*20+30)+10);
}
//---------------------------------------------------------------------------
void __fastcall TForm1::FormKeyDown(TObject *Sender, WORD &Key,
      TShiftState Shift)  {
  switch(Key)  {
  case VK_UP:
    if(dx != 0 && dy != 1)  {
      dx = 0;
      dy = -1;
      }
    break;
  case VK_DOWN:
    if(dx != 0 && dy != -1)  {
      dx = 0;
      dy = 1;
      }
    break;
  case VK_LEFT:
    if(dx != 1 && dy != 0)  {
      dx = -1;
      dy = 0;
      }
    break;
  case VK_RIGHT:
    if(dx != -1 && dy != 0)  {
      dx = 1;
      dy = 0;
      }
      break;
  }
}
//---------------------------------------------------------------------------
void __fastcall TForm1::Timer1Timer(TObject *Sender)  {
  for(int i = snake.size()-1;i >= 0;i--)  {
    if(snake[snake.size()-1].x == apple.x && snake[snake.size()-1].y == apple.y)  {
      snake.push_back(apple);
      apprand();
      count ++;
    }
    if(i != 0)  {
      snake[i].x = snake[i-1].x;
      snake[i].y = snake[i-1].y;
      }
    else  {
      snake[i].x = snake[i].x + dx;
      snake[i].y = snake[i].y + dy;
      }
  }

  for(int i = 1;i < snake.size();i ++)
    if(snake[0].x == snake[i].x && snake[0].y == snake[i].y)  {
      dx = -1;
      dy = count = 0;
      Form1->FormCreate(Form1);
    }
  if(snake[0].x > (x-60)/20 || snake[0].y > (y - 60)/20 || snake[0].x < 0 || snake[0].y < 0)  {
   dx = -1;
   dy = count = 0;
   Form1->FormCreate(Form1);
  }
  Form1->Invalidate();
}
//---------------------------------------------------------------------------

