#pragma once
class Buttonx
{
private:
	int x,y,nomer;
	char* name;
public:
	Buttonx(char*,int,int,int);
	
	void gotoxy(int,int);
	void drav(int);
	void action();
	void SetX(int);
	void SetY(int);
	~Buttonx();
};