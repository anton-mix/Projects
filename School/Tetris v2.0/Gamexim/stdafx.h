// stdafx.h: включаемый файл дл€ стандартных системных включаемых файлов
// или включаемых файлов дл€ конкретного проекта, которые часто используютс€, но
// не часто измен€ютс€
//

#pragma once

#include "targetver.h"
#include <stdio.h>
#include <tchar.h>
#include "conio.h"
#include <stdlib.h>
#include <Windows.h>
#include "Figur.h"
#include <time.h>
#include <ctime>
#include <string.h>
#include <algorithm>
#include "Tetris.h"
#include "Buttonx.h"


// TODO: ”становите здесь ссылки на дополнительные заголовки, требующиес€ дл€ программы
/*#include "stdafx.h"
class Button
{
private:
	int x,y,nomer;
	char* name;
public:
	Button(char* name,int nomer,int x,int y):name(name),nomer(nomer),x(x),y(y){};

	void gotoxy(int xpos, int ypos)	
	{
		COORD scrn;    
		HANDLE hOuput = GetStdHandle(STD_OUTPUT_HANDLE);
		scrn.X = xpos; scrn.Y = ypos;
		SetConsoleCursorPosition(hOuput,scrn);	
	}


	void drav(int nomer1)
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
	
	void action();
		
	~Button(){}
};*/