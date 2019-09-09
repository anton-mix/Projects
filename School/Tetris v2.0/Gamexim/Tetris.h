#pragma once
class Tetris
{
private:
	char pole [22][17],dubpole[20][10],colorpole[22][17];
	
	Figur* activ;
    Figur* passiv;
    Figur* memor;
    float end;
    int score,a;
    bool mem;

	void vvod();
	void clr();
	void update();
	void init();
	void povlev();
	void povprav();
	void povfig();
	bool povniz();
	void memorf();
	void sleep(float);
	void SetColor(char);
	void gotoxy(int,int);
public:
	Tetris();
	~Tetris();
	int Main();
};

