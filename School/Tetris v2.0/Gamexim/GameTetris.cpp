#include "StdAfx.h"
#include "GameTetris.h"

GameTetris::GameTetris()
{
	Game1 = new Tetris();
	cnopca = new Buttonx* [3];
	nomer = 1;
}


GameTetris::~GameTetris()
{
	for(int i = 0;i < 3;i ++)
		delete[] cnopca[i];
	delete [] cnopca;
	delete Game1;

}

void GameTetris::gotoxy(int xpos, int ypos)	
{
	COORD scrn;    
	HANDLE hOuput = GetStdHandle(STD_OUTPUT_HANDLE);
	scrn.X = xpos; scrn.Y = ypos;
	SetConsoleCursorPosition(hOuput,scrn);	
}

void GameTetris::size()
{
	HANDLE hWndConsole;
	if (hWndConsole = GetStdHandle(-12))
    {
        CONSOLE_SCREEN_BUFFER_INFO consoleInfo;
        if (GetConsoleScreenBufferInfo(hWndConsole, &consoleInfo))
        {
            widht = consoleInfo.srWindow.Right - consoleInfo.srWindow.Left + 1;
            height = consoleInfo.srWindow.Bottom - consoleInfo.srWindow.Top + 1;
        }
	}
}

void GameTetris::vvod()
{
	size();
	system("cls");
	
	for(int i = 0;i< 3;i ++)
	{
		cnopca[i]->SetX(widht/3 + 3);
		cnopca[i]->SetY(i * 3 + 2);
	}
	gotoxy((widht)/3 + 1,1);
		
	putchar(218);
	for(int i = 0;i < 14; i ++)
		putchar(196);
	putchar(191);
		
	for(int i = 0;i < 9;i ++)	
	{
		gotoxy((widht)/3 + 1,2 + i);
		putchar(179);
		printf("              ");
		putchar(179);	
	}
		
	gotoxy((widht)/3 + 1,11);
		
	putchar(192);
	for(int i = 0;i < 14; i ++)
		putchar(196);	
	putchar(217);
}

int GameTetris::Main()
{
	size();
	cnopca[0] = new Buttonx("Continue",1,(widht/3 + 3),2);
	cnopca[1] = new Buttonx("Restart",2,(widht/3 + 3),5);
	cnopca[2] = new Buttonx("Exit",3,(widht/3 + 3),8);
	while(true)
	{
	Game1->Main();
	vvod();
	for(int i = 0;i < 3;i ++)
		cnopca[i]->drav(nomer);
l:
	switch(getch())
	{
	case 72:
		nomer = (nomer == 1)?3:nomer - 1;
			for(int i = 0;i < 3;i ++)
				cnopca[i]->drav(nomer);
			goto l;
		break;
	case 80:
		nomer = (nomer == 3)?1:nomer + 1;
			for(int i = 0;i < 3;i ++)
				cnopca[i]->drav(nomer);
			goto l;
		break;
	case 13:
		if(nomer == 1)
		{
			Game1->Main();
			break;
		}
		if(nomer == 2)
		{
			delete Game1;
			Game1 = NULL;
			Game1 = new Tetris();
			break;
		}
		if(nomer == 3)
			return 0;
	case 27:break;
	default:goto l; 
 
  }
	}
	

	return 0;
}