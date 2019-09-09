// Gamexim.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include "GameTetris.h"

int _tmain(int argc, _TCHAR* argv[])
{
	GameTetris* Game = new GameTetris();
	Game->Main();
}

