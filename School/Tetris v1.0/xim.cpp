// xim.cpp: определяет точку входа для консольного приложения.
//


char pole [22][17],dubpole[20][10],colorpole[22][17];
figura* activ = NULL;
figura* passiv = NULL;
figura* memor = NULL;

float end = clock() + 1 * CLOCKS_PER_SEC;
int score = 0;
bool mem = true;

void vvod ()
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

	for(int i = 0;i < 20;i ++)
		for(int j = 0;j < 10;j ++)
			dubpole[i][j] = ' ';
}

void clr()
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

void update()
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
				if(activ->figur[activ->z][i - activ->forynach()][j - activ->forxnach()] == '*')
						pole[i][j] = activ->figur[activ->z][i - activ->forynach()][j - activ->forxnach()];
		}
		//memor
		if(memor != NULL)
		for(int i = memor->forynach();i < memor->forykon();i ++)
			for(int j = memor->forxnach();j < memor->forxkon();j ++)
				pole[i][j] = memor->figur[memor->z][i - memor->forynach()][j - memor->forxnach()];

		//score
		char scor[5];
		sprintf(scor, "%d", score);
		for(int i = 4;i > 0;i --)
		{
			pole[19][i] = scor[4 - i];
			if(scor[4 - i] == '\0')break;
		}
		//pole
		for(int i = 0; i < 22;i ++)
		{
			for(int j  = 0;j < 17;j ++)
				putchar (pole[i][j]);
			putchar('\n');
		}
}

void init()
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
	/*switch (passiv->n)
{
case 1: activ = new line(1,10,0,1);break;
case 10: activ = new cube(1,10,0,10);break;
case 20: activ = new zakprav(1,10,0,20);break;
case 30: activ = new zaklev(1,10,0,30);break;
case 40: activ = new Tfig(1,10,0,40);break;
case 50: activ = new glev(1,10,0,50);break;
case 60: activ = new gprav(1,10,0,60);break;
}*/
std::swap(activ,passiv);
//activ = new figura(passiv);
	activ->pubx(10);
	activ -> puby(1);
	goto l;
}
}

void povlev()
{
	bool prov = true;
for(int i = activ->forynach();i < activ->forykon(); i ++)
{
	for (int j = 0;j < 3; j ++)
		if(activ->figur[activ->z][i - activ->forynach()][j] == '*')
			if(pole[i][activ->forxnach() + (j - 1)] == '*' || activ->forxnach() + (j - 1) == 5)
				prov = false;
			else
				break;
}
if(prov)
	activ->pubx((activ->pubx() - 1));
}

void povprav()
{
	bool prov = true;
for(int i = activ->forynach();i < activ->forykon(); i ++)
{
	for(int j = 0;j > -3;j --)
	if (activ->figur[activ->z][i - activ->forynach()][((activ->forxkon() - 1) - activ->forxnach()) + j] == '*')
		if(pole[i][(activ->forxkon() - 1) + (j + 1)] == '*' || (activ->forxkon() - 1) + (j + 1) == 16)
		prov = false;
		else
			break;
}
if(prov)
	activ->pubx((activ->pubx() + 1));
}

void povfig()
{
	for(int i = activ->forynach();i < activ->forykon();i ++)
		{
			if (i < 1)
				continue;
			for(int j = activ->forxnach();j < activ->forxkon();j ++)
				if(activ->figur[activ->z][i - activ->forynach()][j - activ->forxnach()] == '*')
						pole[i][j] = ' ';
		}

	activ->z = (activ->z == 3)?0:activ->z + 1;
	bool prov = true;
	for(int i = activ->forynach();i < activ->forykon();i ++)
		{
			if (i < 1)
				continue;
			for(int j = activ->forxnach();j < activ->forxkon();j ++)
				if(activ->figur[activ->z][i - activ->forynach()][j - activ->forxnach()] == '*')
						if(pole[i][j] == '*' || pole[i][j] == '+')
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

bool povniz()
{
	bool prov = true;
	for(int i = activ->forxnach();i < activ->forxkon(); i ++)
	{
		for(int j = 0;j > -3;j --)
			if (activ->figur[activ->z][((activ->forykon() - 1) - activ->forynach()) + j][i - activ->forxnach()] == '*')
				if(pole[(activ->forykon() - 1) +(j + 1)][i] == '*' || (activ->forykon() - 1) + (j + 1) == 21 )
					prov = false;
				else
					break;
	}
	return prov;
}

 void memorf()
 {
	 int fign = 0;
	 if (memor != NULL)
	 {
		 fign = memor -> n;
	     delete memor;
	 }

		 switch (activ -> n)
	 {
	 case 1: memor = new line(4,3,0,1);break;
	 case 10: memor = new cube(4,3,0,10);break;
	 case 20: memor = new zakprav(3,3,0,20);break;
	 case 30: memor = new zaklev(3,3,0,30);break;
	 case 40: memor = new Tfig(4,3,0,40);break;
	 case 50: memor = new glev(4,3,0,50);break;
	 case 60: memor = new gprav(4,3,0,60);break;
	 }
		 delete activ;
		 activ = NULL;

	 if (fign == 0)
		  init();
	 else
		 switch (fign)
	 {
	 case 1: activ = new line(1,11,0,1);break;
	 case 10: activ = new cube(1,11,0,10);break;
	 case 20: activ = new zakprav(1,11,0,20);break;
	 case 30: activ = new zaklev(1,11,0,30);break;
	 case 40: activ = new Tfig(1,11,0,40);break;
	 case 50: activ = new glev(1,11,0,50);break;
	 case 60: activ = new gprav(1,11,0,60);break;
	 }

	 mem = false;
 }
void sleep(float seconds = 1)
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
		 {
			 if(dubpole[i][j] == '*')
				 kol ++;
		 }
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

int main()
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

