//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include <stdio.h>
#include "Unit1.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma link "XPMan"
#pragma resource "*.dfm"
TForm1 *Form1;
//---------------------------------------------------------------------------
__fastcall TForm1::TForm1(TComponent* Owner)
        : TForm(Owner)
{
}
//---------------------------------------------------------------------------

void __fastcall TForm1::btn1Click(TObject *Sender)  {
        if(edtInput->Text == '0')
           edtInput->Text = "";
        edtInput->Text = edtInput->Text + ((TButton*)Sender)->Caption;
}
//---------------------------------------------------------------------------
char action;
double n1;
void __fastcall TForm1::btn_pClick(TObject *Sender)  {

   action = ((TButton*)Sender)->Caption.c_str()[0];
   sscanf(edtInput->Text.c_str(),"%lf",&n1);
   edtInput->Text = '0';
}
//---------------------------------------------------------------------------

void __fastcall TForm1::btn_rClick(TObject *Sender)  {
  double n2,result;
  sscanf(edtInput->Text.c_str(),"%lf",&n2);
  switch(action)  {
  case '+':
  result = n1 + n2;
  break;
  case '-':
  result = n1 - n2;
  break;
  case '*':
  result = n1 * n2;
  break;
  case '/':
  if(n2 == 0)  {
   ShowMessage("На 0 делить нельзя");
   return;
  }
  result = n1 / n2;
  break;
  }
   edtInput->Text = result;
}
//---------------------------------------------------------------------------

