#include "StdAfx.h"
#include "Tetris.h"


void Tetris::gotoxy(int xpos, int ypos)	
{
	COORD scrn;    
	HANDLE hOuput = GetStdHandle(STD_OUTPUT_HANDLE);
	scrn.X = xpos; scrn.Y = ypos;
	SetConsoleCursorPosition(hOuput,scrn);	
}


	Tetris::Tetris()
	{
	pole [22][17],dubpole[20][10],colorpole[22][17];
	activ = NULL;
    passiv = NULL;
    memor = NULL;

    end = clock() + 1 * CLOCKS_PER_SEC;
    score = 0;
    mem = true;
	a = 2;
	}

	Tetris::~Tetris()
	{
	delete activ;
	delete passiv;
	delete memor;
	}

	void Tetris::SetColor(char simvol)
	{
		int bagron = 0;
		if(simvol >= 48 && simvol <= 57)
			bagron = 7;
		else
		switch(simvol)
		{
		case '+':bagron = 7;break;
		case 'a':bagron = 14;break;
        case 'b':bagron = 13;break;
        case 'c':bagron = 12;break;
        case 'd':bagron = 11;break;
        case 'e':bagron = 10;break;
        case 'f':bagron = 9;break;
        case 'g':bagron = 15;break;
		}
		HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
		SetConsoleTextAttribute(hConsole,(WORD)(0 << 4)|bagron);
		if(simvol >= 48 && simvol <= 57)
		putch(simvol);
		else
		{
			if(simvol != '+')	
				putch ('*');
			else
				putch ('+');
		}
	}

	void Tetris::vvod()
	{
		for (int i = 0; i < 22;i ++)
			for (int j = 0;j < 17;j ++)
			{
				if(i == 0 ||i == 21)
					pole[i][j] = '+';
				if (j == 0 || j == 16 || j == 5)
					pole[i][j] = '+';
				if((i == 1 || i == 6 || (i >= 11  && i != 19)) && j < 6)
					pole[i][j] = '+';	
			}
			if(passiv == NULL)
			for(int i = 0;i < 20;i ++)
				for(int j = 0;j < 10;j ++)
					dubpole[i][j] = ' ';
	}

	void Tetris::clr()
	{
		system("cls");

		for(int i = 1; i < 21;i ++)
			for(int j  = 6;j < 16;j ++)
				pole[i][j] = dubpole[i - 1][j - 6];
	
		for(int i = 7;i < 11;i ++)
			for(int j = 1;j < 5;j ++)
			{
				pole[i][j] = ' ';
				pole[i - 5][j] = ' ';
				pole[19][j] = ' ';
			}
	}

	void Tetris::povlev()
	{
		bool prov = true;
		for(int i = activ->forynach();i < activ->forykon(); i ++)
		{
			for (int j = 0;j < 3; j ++)
				if(activ->figur[activ->z][i - activ->forynach()][j] == activ->simvol)
				{
					if(pole[i][activ->forxnach() + (j - 1)] != ' ' || activ->forxnach() + (j - 1) == 5)
						prov = false;
					else
						break;
				}
		}
		if(prov)
			activ->pubx((activ->pubx() - 1));
	}


	void Tetris::povprav()
	{
		bool prov = true;
		for(int i = activ->forynach();i < activ->forykon(); i ++)
		{
			for(int j = 0;j > -3;j --)
				if (activ->figur[activ->z][i - activ->forynach()][((activ->forxkon() - 1) - activ->forxnach()) + j] == activ->simvol)
				{
					if(pole[i][(activ->forxkon() - 1) + (j + 1)] != ' ' || (activ->forxkon() - 1) + (j + 1) == 16)
						prov = false;
					else
						break;	
				}
		}
		if(prov)
			activ->pubx((activ->pubx() + 1));
	}


	void Tetris::povfig()
	{
		for(int i = activ->forynach();i < activ->forykon();i ++)
		{
			if (i < 1)	
				continue;
			for(int j = activ->forxnach();j < activ->forxkon();j ++)
				if(activ->figur[activ->z][i - activ->forynach()][j - activ->forxnach()] == activ->simvol)
						pole[i][j] = ' ';
		}
		activ->z = (activ->z == 3)?0:activ->z + 1;
		bool prov = true;
		for(int i = activ->forynach();i < activ->forykon();i ++)
		{
			if (i < 1)
				continue;
			for(int j = activ->forxnach();j < activ->forxkon();j ++)
				if(activ->figur[activ->z][i - activ->forynach()][j - activ->forxnach()] == activ->simvol)
					if(pole[i][j] != ' ' || pole[i][j] == '+')
						prov = false;	
		}

		if(prov == false)
		{
			if (activ->z == 0)
				activ->z = 3;
			else
				activ->z = activ->z - 1;
		}
	}


	bool Tetris::povniz()
	{
		bool prov = true;
		for(int i = activ->forxnach();i < activ->forxkon(); i ++)
		{
			for(int j = 0;j > -3;j --)
				if (activ->figur[activ->z][((activ->forykon() - 1) - activ->forynach()) + j][i - activ->forxnach()] == activ->simvol)
				{
					if(pole[(activ->forykon() - 1) +(j + 1)][i] != ' ' || (activ->forykon() - 1) + (j + 1) == 21 )
						prov = false;
					else
						break;	
				}
		}
		return prov;
	}

	void Tetris::update()
	{
		clr();
		//passiv
		for(int i = passiv->forynach();i < passiv->forykon();i ++)
			for(int j = passiv->forxnach();j < passiv->forxkon();j ++)
				pole[i][j] = passiv->figur[passiv->z][i - passiv->forynach()][j - passiv->forxnach()];	
		//activ
		for(int i = activ->forynach();i < activ->forykon();i ++)
		{
			if (i < 1)
				continue;
			for(int j = activ->forxnach();j < activ->forxkon();j ++)
				if(activ->figur[activ->z][i - activ->forynach()][j - activ->forxnach()] == activ->simvol)
						pole[i][j] = activ->figur[activ->z][i - activ->forynach()][j - activ->forxnach()];
		}
		//memor
		if(memor != NULL)
		for(int i = memor->forynach();i < memor->forykon();i ++)
			for(int j = memor->forxnach();j < memor->forxkon();j ++)
				pole[i][j] = memor->figur[memor->z][i - memor->forynach()][j - memor->forxnach()];

		//score
		char scor[4];
		int a = 4;
		sprintf(scor, "%d", score);
		for(int i = 4;i >= 0;i --)
		{
			if(!(scor[i] >= 48 && scor[i] <= 57))continue;
			pole[19][a] = scor[i];
			a --;
		}
		//pole
		int widht,height;
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
		for(int i = 0; i < 22;i ++)
			{
				gotoxy((widht-1)/3,1 + i);
				for(int j  = 0;j < 17;j ++)	
					SetColor(pole[i][j]);
				putchar('\n');
		}
	}

	void Tetris::init()
	{	
		srand(time(0));
l:
		if(passiv == NULL)
			switch (rand() % 61)
		{
			case 1: passiv = new line(9,3,0,1);break;
			case 10: passiv = new cube(9,3,0,10);break;
			case 20: passiv = new zakprav(8,3,0,20);break;
			case 30: passiv = new zaklev(8,3,0,30);break;
			case 40: passiv = new Tfig(9,3,0,40);break;
			case 50: passiv = new glev(9,3,0,50);break;
			case 60: passiv = new gprav(9,3,0,60);break;
			default: goto l;
		}
		if(activ == NULL)
		{
			std::swap(passiv,activ);
			activ->pubx(10);
			activ->puby(1);
			delete passiv;
			passiv = NULL;
			goto l;
		}
	}

	void Tetris::memorf()
 {
	 Figur* buff = NULL;
	 std::swap(memor,buff);
	 std::swap(activ,memor);
	 std::swap(buff,activ);
	 if(activ == NULL)
		 init();
	 activ->pubx(10);
	 activ->puby(1);
	 memor->pubx(3);
	 memor->puby(4);

	 mem = false;
	 update();
 }

	void Tetris::sleep(float seconds = 1)
	{
		if (clock() > end)
		{
			if(povniz())
				activ->puby(activ->puby() + 1);
			else
			{
				for(int i = 1; i < 21;i ++)
					for(int j  = 6;j < 16;j ++)
						dubpole[i - 1][j - 6] = pole[i][j];

				mem = true;
				delete activ;
				activ = NULL;
				init();
				
				for(int i = 0;i < 20;i ++)
				{
					int kol = 0;
					for(int j = 0;j < 10;j ++)
						if(dubpole[i][j] != ' ')
							kol ++;
					
					if(kol == 10)
					{
						for(int l = 0;l < 10;l ++)
							dubpole[i][l] = ' ';
						for(int y = i - 1;y > 0;y --)
							for(int x = 0;x < 10 ;x ++)
								dubpole[y + 1][x] = dubpole[y][x];
						i --;
						score ++;
					}
				}
			}

			end = clock() + seconds * CLOCKS_PER_SEC;
		}
	}

	int Tetris::Main()
	{
		vvod();
		init();
		sleep();
		update();
		int x = activ->pubx();
		int y = activ->puby();
		int z = activ->z;
		while(true)
		{
			a = a;
		if(kbhit())
		{
			switch(getch())
			{
			case 75:povlev();break;
			case 77:povprav();break;
			case 72:
				if(mem)
					memorf();
				break;
			case 80:
				if(povniz())
					activ->puby(activ->puby() + 1);
				break;
			case 13:povfig();break;
			case 27:return 0;
			}
		}
		

		sleep();

		if(activ->pubx() != x || activ->puby() != y || activ->z != z)
		{
			update();
			x = activ->pubx();
			y = activ->puby();
			z = activ->z;
		
		}
		}
		return 0;
}






