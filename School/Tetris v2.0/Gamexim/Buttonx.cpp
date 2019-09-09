#include "StdAfx.h"
#include "Buttonx.h"

Buttonx::Buttonx(char* name,int nomer,int x,int y):name(name),nomer(nomer),x(x),y(y)
{}

Buttonx::~Buttonx(void)
{}

void Buttonx::gotoxy(int xpos, int ypos)	
	{
		COORD scrn;    
		HANDLE hOuput = GetStdHandle(STD_OUTPUT_HANDLE);
		scrn.X = xpos; scrn.Y = ypos;
		SetConsoleCursorPosition(hOuput,scrn);	
	}

void Buttonx::SetX(int x)
{
	this->x = x;
}

void Buttonx::SetY(int y)
{
	this->y = y;
}

void Buttonx::drav(int nomer1)
	{
		if(nomer1 != nomer)
		{
			gotoxy(x,y);
		
			putchar(218);
			for(int i = 0;i< 10;i ++)
				putchar(196);
			putchar(191);
		
		
			gotoxy(x,y+1);
		
			putchar(179);
			printf("%s",name);
			for(int i = 0;i < 10 - strlen(name);i ++)
				putchar(' ');
			putchar(179);
		
			gotoxy(x,y+2);
		
			putchar(192);
			for(int i = 0;i < 10;i ++)
				putchar(196);
			putchar(217);
		}
		else
		{
		
			gotoxy(x,y);
		
			putchar(201);
			for(int i = 0;i< 10;i ++)
				putchar(205);
			putchar(187);
		
		
			gotoxy(x,y+1);
		
			putchar(186);
			printf("%s",name);
			for(int i = 0;i < 10 - strlen(name);i ++)
				putchar(' ');
			putchar(186);
		
			gotoxy(x,y+2);
		
			putchar(200);
			for(int i = 0;i < 10;i ++)
				putchar(205);
			putchar(188);
		}
	}

