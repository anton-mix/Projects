//---------------------------------------------------------------------------

#ifndef Unit1H
#define Unit1H
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <ExtCtrls.hpp>
#include "XPMan.hpp"
//---------------------------------------------------------------------------
class TForm1 : public TForm
{
__published:	// IDE-managed Components
        TButton *btn1;
        TPanel *Panel1;
        TButton *btn2;
        TButton *btn3;
        TButton *btn4;
        TButton *btn5;
        TButton *btn6;
        TEdit *edtInput;
        TButton *btn7;
        TButton *btn8;
        TButton *btn9;
        TButton *btn_t;
        TButton *btn0;
        TButton *btn_r;
        TButton *btn_p;
        TButton *btn_m;
        TButton *btn_d;
        TButton *btn_u;
        void __fastcall btn1Click(TObject *Sender);
        void __fastcall btn_pClick(TObject *Sender);
        void __fastcall btn_rClick(TObject *Sender);
private:	// User declarations
public:		// User declarations
        __fastcall TForm1(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TForm1 *Form1;
//---------------------------------------------------------------------------
#endif
