#pragma once
class GameTetris
{
	int widht,height;
	Tetris* Game1;
    Buttonx** cnopca;
	int nomer;
	void gotoxy(int,int);	
	void size();
	void vvod();
public:
	GameTetris();
	~GameTetris();
	int Main();
};

