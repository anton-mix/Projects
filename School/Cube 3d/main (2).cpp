#include <Windows.h>
#include <math.h>
#include <vector>

struct vertex
{
double x,y,z;
vertex(double x = 0,double y = 0,double z = 0):x(x),y(y),z(z){};
};
std::vector<vertex> m_vertex;
typedef std::vector<int> vect;
std::vector<vect> m_faces;

BOOL InitApp(HINSTANCE);
LRESULT CALLBACK WndProc(HWND, UINT, WPARAM, LPARAM);

WCHAR const szClassName[] = L"WindowAppClass";

WCHAR const szWindowTitle[] = L"asfd";

int WinMain(HINSTANCE hInstance,HINSTANCE hPrevInstance,LPSTR lpszCmdLine,int nCmdShow)
{
MSG msg;
HWND hwnd;
if(!hPrevInstance)
{
	if (!InitApp(hInstance))
		return FALSE;
}
else
{
	MessageBox(NULL,(LPCSTR)"Можно запускать только одну копию приложения",(LPCSTR)"Ошибка",MB_OK | MB_ICONSTOP);
	return FALSE;
}

hwnd = CreateWindow((LPCSTR)szClassName,(LPCSTR)szWindowTitle,WS_MINIMIZEBOX|WS_CAPTION|WS_SYSMENU,500,200,700,500,0,0,hInstance,NULL);

if(!hwnd)
	return FALSE;

ShowWindow(hwnd,nCmdShow);
UpdateWindow(hwnd);

while (GetMessage(&msg, 0, 0, 0))
DispatchMessage(&msg);

return msg.wParam;

}

BOOL InitApp(HINSTANCE hInstance)
{
ATOM aWndClass;
WNDCLASS wc;

memset (&wc, 0, sizeof(wc));
wc.style = CS_HREDRAW | CS_VREDRAW;

wc.lpfnWndProc = (WNDPROC)  WndProc;

wc.cbClsExtra = 0;

wc.cbWndExtra = 0;

wc.hInstance = hInstance;

wc.hIcon = LoadIcon(NULL, IDI_APPLICATION);

wc.hCursor = LoadCursor(NULL, IDC_ARROW);

wc.hbrBackground = CreateSolidBrush(RGB(0x67,0x8F,0xC2));

wc.lpszMenuName = (LPCSTR)NULL;

wc.lpszClassName = (LPCSTR)szClassName;

aWndClass = RegisterClass(&wc);

return (aWndClass != 0);
}

void fillvertex(int width)
{
m_vertex.resize(8);
m_vertex[0] = vertex(-width/2,-width/2,width/2);
m_vertex[1] = vertex(width/2,-width/2,width/2);
m_vertex[2] = vertex(width/2,width/2,width/2);
m_vertex[3] = vertex(-width/2,width/2,width/2);
m_vertex[4] = vertex(-width/2,-width/2,-width/2);
m_vertex[5] = vertex(width/2,-width/2,-width/2);
m_vertex[6] = vertex(width/2,width/2,-width/2);
m_vertex[7] = vertex(-width/2,width/2,-width/2);
}

void fillfaces()
{
m_faces.resize(6);
for(int i = 0;i < 6;i ++)
	m_faces[i].resize(4);
m_faces[0][0] = 0;m_faces[0][1] = 1;m_faces[0][2] = 2;m_faces[0][3] = 3;
m_faces[1][0] = 4;m_faces[1][1] = 5;m_faces[1][2] = 6;m_faces[1][3] = 7;
m_faces[2][0] = 0;m_faces[2][1] = 1;m_faces[2][2] = 5;m_faces[2][3] = 4;
m_faces[3][0] = 1;m_faces[3][1] = 2;m_faces[3][2] = 6;m_faces[3][3] = 5;
m_faces[4][0] = 2;m_faces[4][1] = 3;m_faces[4][2] = 7;m_faces[4][3] = 6;
m_faces[5][0] = 3;m_faces[5][1] = 0;m_faces[5][2] = 4;m_faces[5][3] = 7;
}

void rotate(double PhiX,double PhiY,double PhiZ)
{
	double ca = cos(PhiX);
	double sa = sin(PhiX);
	double cb = cos(PhiY);
	double sb = sin(PhiY);
	double cg = cos(PhiZ);
	double sg = sin(PhiZ);

	for(int i = 0; i < (int)m_vertex.size();i ++)
	{
	double x = m_vertex[i].x;
	double y = m_vertex[i].y;
	double z = m_vertex[i].z;
	double x2 = (cb*cg) * x - (cb*sg) * y + sb * z;
	double y2 = (cg*sa*sb + ca*sg) * x + (ca*cg - sa*sb*sg) * y - (cb*sa) * z;
	double z2 = (-ca*cg*sb + sa*sg) * x + (cg*sa + ca*sb*sg) * y + (ca*cb) * z;
	m_vertex[i] = vertex(x2,y2,z2);
	}
}

void draw(HDC hdc,int x0,int y0)
{
POINT pts[4];
for(int i = 0; i < (int)m_faces.size();i ++)
{
	for(int j = 0;j < m_faces[i].size();j ++)
	{
		pts[j].x = m_vertex[m_faces[i][j]].x + x0;
		pts[j].y = y0 - m_vertex[m_faces[i][j]].y;
	}
	Polygon(hdc,pts,m_faces[i].size());
}
}

LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam,LPARAM lParam)
{
PAINTSTRUCT ps;
HDC hdc;
RECT rc;
GetClientRect(hWnd,&rc);
switch(message)
{
case WM_ACTIVATE:
	fillvertex(100);
	fillfaces();
	break;
case WM_PAINT:
	hdc = BeginPaint(hWnd, &ps);
	draw(hdc,(rc.right - rc.left)/2,(rc.bottom - rc.top) / 2);
	EndPaint(hWnd, &ps);
	break;
case WM_KEYDOWN:
	switch(wParam)
	{
	case VK_LEFT:
		rotate(0,0,-5.0 * 3.14/180.0);
		break;
	case VK_RIGHT:
		rotate(0,0,5.0 * 3.14/180.0);
		break;
	case VK_UP:
		rotate(0,5.0*3.14/180.0,0);
		break;
	case VK_DOWN:
		rotate(0,-5.0*3.14/180.0,0);
		break;
	}
    InvalidateRect(hWnd,NULL,true);
    break;
case WM_SIZE:
	InvalidateRect(hWnd,NULL,true);
	break;
case WM_DESTROY:
	PostQuitMessage(0);
	break;
default:
		return DefWindowProc(hWnd,message,wParam,lParam);
}
return 0;
}
