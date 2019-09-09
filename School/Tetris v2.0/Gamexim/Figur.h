#pragma once
class Figur
{
protected:
	int x,y;

public:
	char figur[4][4][4],simvol;
	int z,n;
	Figur(int y1,int x1,int z1,int n1)
	{
		x = x1;
		y = y1;
		z = z1;
		n = n1;
	}
	
	Figur(){}
	
	virtual int forynach() = 0;
	virtual int forykon() = 0;
	virtual int forxnach() = 0;
	virtual int forxkon() = 0;


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

class line : public Figur
{
public:
	
	line (int y1,int x1,int z1,int n1) : Figur(y1,x1,z1,n1)
	{
		simvol = 'a';
		for(int i = 0;i < 4;i ++)
		{
			figur[0][i][0] = simvol;
			figur[1][0][i] = simvol;
			figur[2][i][0] = simvol;
			figur[3][0][i] = simvol;
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

class cube : public Figur
{
public :
	cube (int y1,int x1,int z1,int n1) : Figur(y1,x1,z1,n1){
	
		simvol = 'b';

		for(int j = 0;j < 4; j ++)
			for(int i = 0;i < 2; i ++)
			{
				figur[j][0][i] = simvol;
				figur[j][1][i] = simvol;
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

class zakprav : public Figur
{
public :
	zakprav (int y1,int x1,int z1,int n1) : Figur(y1,x1,z1,n1)
	{
		simvol = 'c';

		for(int i = 0;i < 3;i ++)
		{
			figur[0][0][i] = (i == 0)?' ':simvol;	
			figur[0][1][i] = (i == 2)?' ':simvol;
	
			figur[1][i][1] = (i == 0)?' ':simvol;
			figur[1][i][0] = (i == 2)?' ':simvol;

			figur[2][0][i] = (i == 0)?' ':simvol;
			figur[2][1][i] = (i == 2)?' ':simvol;

			figur[3][i][1] = (i == 0)?' ':simvol;
			figur[3][i][0] = (i == 2)?' ':simvol;
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

class zaklev : public Figur
{
public :
	zaklev (int y1,int x1,int z1,int n1) : Figur(y1,x1,z1,n1)
	{
		simvol = 'd';

		for(int i = 0;i < 3;i ++)
		{
	
			figur[0][0][i] = (i == 2)?' ':simvol;
			figur[0][1][i] = (i == 0)?' ':simvol;
	
			figur[1][i][1] = (i == 2)?' ':simvol;
			figur[1][i][0] = (i == 0)?' ':simvol;
	
			figur[2][0][i] = (i == 2)?' ':simvol;
			figur[2][1][i] = (i == 0)?' ':simvol;
	
			figur[3][i][1] = (i == 2)?' ':simvol;
			figur[3][i][0] = (i == 0)?' ':simvol;
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

class Tfig : public Figur
{
public :
	Tfig (int y1,int x1,int z1,int n1) : Figur(y1,x1,z1,n1)
	{
		simvol = 'e';

		for(int i = 0;i < 3;i ++)
		{
			figur[0][0][i] = (i == 1)?simvol:' ';
			figur[0][1][i] = simvol;
	
			figur[1][i][1] = (i == 1)?simvol:' ';
			figur[1][i][0] = simvol;
	
			figur[2][1][i] = (i == 1)?simvol:' ';
			figur[2][0][i] = simvol;

			figur[3][i][0] = (i == 1)?simvol:' ';
			figur[3][i][1] = simvol;
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

class glev : public Figur
{
public :
	glev (int y1,int x1,int z1,int n1) : Figur(y1,x1,z1,n1)
	{
		simvol = 'f';
		
		for(int i = 0;i < 3;i ++)
		{
		figur[0][1][i] = simvol;
		figur[0][0][i] = (i == 2)?simvol:' ';
	
	
		figur[1][i][0] = simvol;
		figur[1][i][1] = (i == 2)?simvol:' ';
	
		figur[2][0][i] = simvol;
		figur[2][1][i] = (i == 0)?simvol:' ';
	
		figur[3][i][1] = simvol;
		figur[3][i][0] = (i == 0)?simvol:' ';
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

class gprav : public Figur
{
public :
	gprav (int y1,int x1,int z1,int n1) : Figur(y1,x1,z1,n1)
	{
		simvol = 'g';

		for(int i = 0;i < 3;i ++)
		{
			figur[0][0][i] = (i == 0)?simvol:' ';
			figur[0][1][i] = simvol;
	
			figur[1][i][1] = (i == 0)?simvol:' ';
			figur[1][i][0] = simvol;

			figur[2][1][i] = (i == 2)?simvol:' ';
			figur[2][0][i] = simvol;

			figur[3][i][0] = (i == 2)?simvol:' ';
			figur[3][i][1] = simvol;
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

