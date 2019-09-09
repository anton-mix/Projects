class figura{
protected:
	int x,y;

public:
	char figur[4][3][3];
	int z,n;
	figura(int y1,int x1,int z1,int n1)
{
	x = x1;
	y = y1;
	z = z1;
	n = n1;
}
figura (const figura& a){
x=a.x;
y=a.y;
z=a.z;
n=a.n;
}
	figura(){}

	virtual int forynach(){};
	virtual int forykon(){};
	virtual int forxnach(){};
	virtual int forxkon(){};

int pubx(int x1 = NULL)
{
	if (x1 != NULL)
		x = x1;
	else
		return x;
}
int puby(int y1 = NULL)
{
	if(y1 != NULL)
		y = y1;
	else
		return y;
}
};

class line : public figura
{
public:
	line (int y1,int x1,int z1,int n1) : figura(y1,x1,z1,n1){
		for(int i = 0;i < 4;i ++)
		{
			figur[0][i][0] = '*';
			figur[1][0][i] = '*';
			figur[2][i][0] = '*';
			figur[3][0][i] = '*';
		}
	}

	int forynach()
	{
		if(z == 0 || z == 2)
	    return y - 2;
		else
		return y;
	}

	int forykon()
	{
		if (z == 0 || z == 2)
		return y + 2;
		else
		return y + 1;
	}
	int forxnach()
	{
		if (z == 0 || z == 2)
		return x;
		else
		return x - 2;
	}

	int forxkon()
	{
		if (z == 0 || z == 2)
		return x + 1;
		else
		return x + 2;
	}
};

class cube : public figura
{
public :
	cube (int y1,int x1,int z1,int n1) : figura(y1,x1,z1,n1){

		for(int j = 0;j < 4; j ++)
			for(int i = 0;i < 2; i ++)
			{
				figur[j][0][i] = '*';
				figur[j][1][i] = '*';
			}
	}

	int forynach()
	{
		return y - 1;
	}

	int forykon()
	{
		return y + 1;
	}

	int forxnach()
	{
		return x - 1;
	}

	int forxkon()
	{
		return x + 1;
	}

};

class zakprav : public figura
{
public :
	zakprav (int y1,int x1,int z1,int n1) : figura(y1,x1,z1,n1){
	for(int i = 0;i < 3;i ++)
	{
	figur[0][0][i] = (i == 0)?' ':'*';
	figur[0][1][i] = (i == 2)?' ':'*';

	figur[1][i][1] = (i == 0)?' ':'*';
	figur[1][i][0] = (i == 2)?' ':'*';

	figur[2][0][i] = (i == 0)?' ':'*';
	figur[2][1][i] = (i == 2)?' ':'*';

	figur[3][i][1] = (i == 0)?' ':'*';
	figur[3][i][0] = (i == 2)?' ':'*';
	}
	}

	int forynach()
	{
		if (z == 0 || z == 2)
			return y;
		else
			return y - 1;

	}

	int forykon()
	{
			return y + 2;
	}

	int forxnach()
	{
		if (z == 0 || z == 2)
			return x - 1;
		else
			return x;
	}

	int forxkon()
	{
			return x + 2;
	}
};

class zaklev : public figura
{
public :
	zaklev (int y1,int x1,int z1,int n1) : figura(y1,x1,z1,n1){
	for(int i = 0;i < 3;i ++)
	{
	figur[0][0][i] = (i == 2)?' ':'*';
	figur[0][1][i] = (i == 0)?' ':'*';

	figur[1][i][1] = (i == 2)?' ':'*';
	figur[1][i][0] = (i == 0)?' ':'*';

	figur[2][0][i] = (i == 2)?' ':'*';
	figur[2][1][i] = (i == 0)?' ':'*';

	figur[3][i][1] = (i == 2)?' ':'*';
	figur[3][i][0] = (i == 0)?' ':'*';
	}
	}

	int forynach()
	{
	if (z == 0 || z == 2)
		return y;
	else
		return y - 1;
	}

	int forykon()
	{
		return y + 2;
	}

	int forxnach()
	{
		return x - 1;
	}

	int forxkon()
	{
		if (z == 0 || z == 2)
			return x + 2;
		else
			return x + 1;
	}
};

class Tfig : public figura
{
public :
	Tfig (int y1,int x1,int z1,int n1) : figura(y1,x1,z1,n1){
	for(int i = 0;i < 3;i ++)
	{
	figur[0][0][i] = (i == 1)?'*':' ';
	figur[0][1][i] = '*';

	figur[1][i][1] = (i == 1)?'*':' ';
	figur[1][i][0] = '*';

	figur[2][1][i] = (i == 1)?'*':' ';
	figur[2][0][i] = '*';

	figur[3][i][0] = (i == 1)?'*':' ';
	figur[3][i][1] = '*';
	}
	}

	int forynach()
	{
		if(z == 2)
			return y;
		else
		    return y - 1;
	}

	int forykon()
	{
		if (z == 0)
			return y + 1;
		else
			return y + 2;
	}

	int forxnach()
	{
		if (z == 1)
			return x;
		else
			return x - 1;
	}

	int forxkon()
	{
		if (z == 3)
			return x + 1;
		else
			return x + 2;
	}
};

class glev : public figura
{
public :
	glev (int y1,int x1,int z1,int n1) : figura(y1,x1,z1,n1){
	for(int i = 0;i < 3;i ++)
	{
	figur[0][1][i] = '*';
	figur[0][0][i] = (i == 2)?'*':' ';

	figur[1][i][0] = '*';
	figur[1][i][1] = (i == 2)?'*':' ';

	figur[2][0][i] = '*';
	figur[2][1][i] = (i == 0)?'*':' ';

	figur[3][i][1] = '*';
	figur[3][i][0] = (i == 0)?'*':' ';
	}
	}

	int forynach()
	{
		if (z == 2)
			return y;
		else
			return y - 1;
	}

	int forykon()
	{
		if(z == 0)
			return y + 1;
		else
			return y + 2;
	}

	int forxnach()
	{
		if(z == 1)
			return x;
		else
			return x - 1;
	}

	int forxkon()
	{
		if(z == 3)
			return x + 1;
		else
			return x + 2;
	}
};

class gprav : public figura
{
public :
	gprav (int y1,int x1,int z1,int n1) : figura(y1,x1,z1,n1){
	for(int i = 0;i < 3;i ++)
	{
	figur[0][0][i] = (i == 0)?'*':' ';
	figur[0][1][i] = '*';

	figur[1][i][1] = (i == 0)?'*':' ';
	figur[1][i][0] = '*';

	figur[2][1][i] = (i == 2)?'*':' ';
	figur[2][0][i] = '*';

	figur[3][i][0] = (i == 2)?'*':' ';
	figur[3][i][1] = '*';
	}
	}

	int forynach()
	{
		if(z == 2)
			return y;
		else
			return y - 1;
	}

	int forykon()
	{
		if(z == 0)
			return y + 1;
		else
			return y + 2;
	}

	int forxnach()
	{
		if(z == 1)
			return x;
		else
			return x - 1;
	}

	int forxkon()
	{
		if(z == 3)
			return x + 1;
		else
			return x + 2;
	}
};
