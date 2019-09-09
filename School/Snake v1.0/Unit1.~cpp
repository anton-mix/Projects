//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "Unit1.h"
#include <time.h>
#include <stdlib.h>
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TForm1 *Form1;
//---------------------------------------------------------------------------
__fastcall TForm1::TForm1(TComponent* Owner)
        : TForm(Owner)
{
}
//---------------------------------------------------------------------------
int width,height,size = 20;
TPoint food;
bool autop;
void size1()  {
 width = Form1->ClientWidth / size;
 height = Form1->ClientHeight / size;
}

void TForm1::genFood()  {
  while(true)  {
  int x = rand()%width;
  int y = rand()%height;
  bool exsis = false;
  for(int i = 0; i < Form1->snake->gettail().size();i ++)  {
   if(exsis = (x == Form1->snake->gettail()[i].x && y == Form1->snake->gettail()[i].y))
     break;
  }
   if(!exsis)  {
   food.x  = x;
   food.y  = y;
   Form1->Invalidate();
   break;
   }
}
}

void __fastcall TForm1::FormCreate(TObject *Sender)
{
  size1();
  Form1->snake = new Snake(width/2,height/2);
  Form1->DoubleBuffered = true;
  srand(time(0));
  Form1->genFood();
}
//---------------------------------------------------------------------------

void __fastcall TForm1::FormResize(TObject *Sender)
{
  size1();
  Form1->Invalidate();
}
//---------------------------------------------------------------------------
void __fastcall TForm1::FormPaint(TObject *Sender)  {
  Form1->Canvas->Brush->Color = clRed;
  Form1->Canvas->Rectangle(0,0,Form1->ClientWidth,Form1->ClientHeight);
  if(Form1->Timer1->Enabled) {
   Form1->Canvas->Brush->Style = bsClear;
   Form1->Canvas->Font->Name = "Arial";
   Form1->Canvas->Font->Size = 40;
   Form1->Canvas->Font->Color = 0xc0c0c0;
   AnsiString count1 = Form1->snake->gettail().size();
   Form1->Canvas->TextOutA(
   (Form1->ClientWidth - Form1->Canvas->TextWidth(count1))/2,
   (Form1->ClientHeight - Form1->Canvas->TextHeight(count1))/2,
    count1);
  }
  int offsetX = (Form1->ClientWidth - width * size) / 2;
  int offsetY = (Form1->ClientHeight - height * size) / 2;
  Form1->Canvas->Pen->Style = psSolid;
  Form1->Canvas->Pen->Color = clBlue;

  for(int i = offsetX;i <= Form1->ClientWidth - offsetX;i += size)  {
   Form1->Canvas->MoveTo(i,offsetY);
   Form1->Canvas->LineTo(i,Form1->ClientHeight - offsetY);
  }

  for(int i = offsetY;i <= Form1->ClientHeight - offsetY; i += size)  {
   Form1->Canvas->MoveTo(offsetX, i);
   Form1->Canvas->LineTo(Form1->ClientWidth - offsetX, i);
  }

  Form1->Canvas->Pen->Style = psClear;
  Form1->Canvas->Brush->Color = clYellow;
  Form1->Canvas->FillRect(TRect(offsetX  + food.x * size + 1,
                                 offsetY + food.y * size + 1,
                                 offsetX + (food.x+1) *size,
                                 offsetY + (food.y+1)*size));

  for(int i = 0;i <Form1->snake->gettail().size();i++)  {
   int x = Form1->snake->gettail()[i].x;
   int y = Form1->snake->gettail()[i].y;
   Form1->Canvas->Pen->Style = psClear;
   Form1->Canvas->Brush->Color = clGreen;
   Form1->Canvas->FillRect(Rect(offsetX + x *size + 1,offsetY + y * size + 1,
   offsetX + (x+1)*size,offsetY+ (y+1) * size));

   if(Form1->Timer1->Enabled == false)  {
   Form1->Canvas->Brush->Style = bsClear;
   Form1->Canvas->Font->Name = "Arial";
   Form1->Canvas->Font->Size = 20;
   Form1->Canvas->Font->Color = 0xc0c0c0;
   AnsiString count = "Game over";
   Form1->Canvas->TextOutA(
   (Form1->ClientWidth - Form1->Canvas->TextWidth(count))/2,
   (Form1->ClientHeight - Form1->Canvas->TextHeight(count))/2,
    count);
   }
  }
}
//---------------------------------------------------------------------------
void __fastcall TForm1::Timer1Timer(TObject *Sender)  {
  if(autop)  {
   if(Form1->snake->gethead().y != food.y)
     Form1->snake->setDerection((food.y > Form1->snake->gethead().y)?bottom:top);
    if(Form1->snake->gethead().x != food.x)
     Form1->snake->setDerection((food.x > Form1->snake->gethead().x)?right:left);
  }
  Form1->Timer1->Enabled = Form1->snake->step();
  if(Form1->snake->gethead().x == food.x && Form1->snake->gethead().y == food.y)  {
  Form1->snake->eat(food);
  Form1->genFood();
  }
  if(Form1->snake->gethead().x < 0  || Form1->snake->gethead().y < 0  ||
     Form1->snake->gethead().x >= width || Form1->snake->gethead().y >= height)  {
   Form1->Timer1->Enabled = false;

   }
  Form1->Invalidate();
}
//---------------------------------------------------------------------------

void __fastcall TForm1::FormKeyDown(TObject *Sender, WORD &Key,
      TShiftState Shift)
{
 switch(Key)  {
  case VK_LEFT:
  if(Form1->snake->getDerection() != right) {
  Form1->snake->setDerection(left);
  }
  break;
  case VK_RIGHT:
  if(Form1->snake->getDerection() != left) {
  Form1->snake->setDerection(right);
  }
  break;
  case VK_UP:
  if(Form1->snake->getDerection() != bottom) {
  Form1->snake->setDerection(top);
  }
  break;
  case VK_DOWN:
  if(Form1->snake->getDerection() != top) {
  Form1->snake->setDerection(bottom);
  }
  break;
  case VK_ESCAPE:
  autop = !autop;
  break;
  }
}
//---------------------------------------------------------------------------

void __fastcall TForm1::FormClick(TObject *Sender)
{
if(!Form1->Timer1->Enabled)  {
 Form1->snake->reset(width/2,height/2);
 Form1->genFood();
 Form1->Timer1->Enabled = true;
 Form1->Invalidate();
}
}
//---------------------------------------------------------------------------

