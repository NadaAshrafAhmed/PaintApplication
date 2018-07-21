// Project.cpp : Defines the entry point for the application.
//
#pragma warning(disable:4996)
#include "stdafx.h"
#include "Project.h"
#include<algorithm>
#include<iostream>
#include <string>
#include <vector>
#define _USE_MATH_DEFINES
#include "math.h"
using namespace std;

#define MAX_LOADSTRING 100

// Global Variables:
HINSTANCE hInst;                                // current instance
WCHAR szTitle[MAX_LOADSTRING];                  // The title bar text
WCHAR szWindowClass[MAX_LOADSTRING];            // the main window class name

// Forward declarations of functions included in this code module:
ATOM                MyRegisterClass(HINSTANCE hInstance);
BOOL                InitInstance(HINSTANCE, int);
LRESULT CALLBACK    WndProc(HWND, UINT, WPARAM, LPARAM);
INT_PTR CALLBACK    About(HWND, UINT, WPARAM, LPARAM);

int APIENTRY wWinMain(_In_ HINSTANCE hInstance,
                     _In_opt_ HINSTANCE hPrevInstance,
                     _In_ LPWSTR    lpCmdLine,
                     _In_ int       nCmdShow)
{
    UNREFERENCED_PARAMETER(hPrevInstance);
    UNREFERENCED_PARAMETER(lpCmdLine);

    // TODO: Place code here.

    // Initialize global strings
    LoadStringW(hInstance, IDS_APP_TITLE, szTitle, MAX_LOADSTRING);
    LoadStringW(hInstance, IDC_PROJECT, szWindowClass, MAX_LOADSTRING);
    MyRegisterClass(hInstance);

    // Perform application initialization:
    if (!InitInstance (hInstance, nCmdShow))
    {
        return FALSE;
    }

    HACCEL hAccelTable = LoadAccelerators(hInstance, MAKEINTRESOURCE(IDC_PROJECT));

    MSG msg;

    // Main message loop:
    while (GetMessage(&msg, nullptr, 0, 0))
    {
        if (!TranslateAccelerator(msg.hwnd, hAccelTable, &msg))
        {
            TranslateMessage(&msg);
            DispatchMessage(&msg);
        }
    }

    return (int) msg.wParam;
}



//
//  FUNCTION: MyRegisterClass()
//
//  PURPOSE: Registers the window class.
//
ATOM MyRegisterClass(HINSTANCE hInstance)
{
    WNDCLASSEXW wcex;

    wcex.cbSize = sizeof(WNDCLASSEX);

    wcex.style          = CS_HREDRAW | CS_VREDRAW;
    wcex.lpfnWndProc    = WndProc;
    wcex.cbClsExtra     = 0;
    wcex.cbWndExtra     = 0;
    wcex.hInstance      = hInstance;
    wcex.hIcon          = LoadIcon(hInstance, MAKEINTRESOURCE(IDI_PROJECT));
    wcex.hCursor        = LoadCursor(nullptr, IDC_ARROW);
    wcex.hbrBackground  = (HBRUSH)(COLOR_WINDOW+1);
    wcex.lpszMenuName   = MAKEINTRESOURCEW(IDC_PROJECT);
    wcex.lpszClassName  = szWindowClass;
    wcex.hIconSm        = LoadIcon(wcex.hInstance, MAKEINTRESOURCE(IDI_SMALL));

    return RegisterClassExW(&wcex);
}

//
//   FUNCTION: InitInstance(HINSTANCE, int)
//
//   PURPOSE: Saves instance handle and creates main window
//
//   COMMENTS:
//
//        In this function, we save the instance handle in a global variable and
//        create and display the main program window.
//
BOOL InitInstance(HINSTANCE hInstance, int nCmdShow)
{
   hInst = hInstance; // Store instance handle in our global variable

   HWND hWnd = CreateWindowW(szWindowClass, szTitle, WS_OVERLAPPEDWINDOW,
      CW_USEDEFAULT, 0, CW_USEDEFAULT, 0, nullptr, nullptr, hInstance, nullptr);

   if (!hWnd)
   {
      return FALSE;
   }

   ShowWindow(hWnd, nCmdShow);
   UpdateWindow(hWnd);

   return TRUE;
}

//Lines

//Parametric

void ParametricLine(HDC hdc, int xs, int ys, int xe, int ye, COLORREF color) {

	int deltax = xe - xs;
	int deltay = ye - ys;

	int n = max(abs(deltax), abs(deltay));

	double dt = 1.0 / n;
	double dx = dt * (double)deltax;
	double dy = dt * (double)deltay;

	double x = xs;
	double y = ys;

	for (int i = 0; i < n; i++) {

		SetPixel(hdc, round(x), round(y), color);

		x += dx;
		y += dy;
	}
}

//DDA

void DDA(HDC hdc, int xs, int ys, int xe, int ye, COLORREF color)
{
	int delta_x = (xe - xs);
	int delta_y = (ye - ys);
	double m = delta_y / (double)delta_x;
	if (abs(delta_y) < abs(delta_x))
	{
		if (xe < xs)
		{
			swap(xs, xe);
			swap(ys, ye);
		}
		int x = xs;
		double y = ys;
		while (x < xe)
		{
			SetPixel(hdc, x, round(y), color);
			x++;
			y += m;
		}
	}
	else
	{
		m = 1.0 / m;
		if (ye < ys)
		{
			swap(xs, xe);
			swap(ys, ye);
		}
		double x = xs;
		int y = ys;
		while (y < ye)
		{
			SetPixel(hdc, round(x), y, color);
			x += m;
			y++;
		}
	}
}

//Midpoint

void LineMidPoint(HDC hdc, int xs, int ys, int xe, int ye, COLORREF color) {

	double slope = (double)(ye - ys) / (xe - xs);

	if (slope < 0) {

		slope = slope *-1;
	}

	if (slope <= 1) {

		if (xs > xe) {

			int temp = xs;
			xs = xe;
			xe = temp;

			temp = ys;
			ys = ye;
			ye = temp;
		}

		int incy;

		if (ys > ye) {
			incy = -1;
		}
		else {
			incy = 1;
		}

		int deltax = xe - xs;
		int deltay = ye - ys;

		int d = deltax - (2 * abs(deltay));

		int d1 = (2 * deltax) - (2 * abs(deltay));
		int d2 = -2 * abs(deltay);

		int x = xs;
		int y = ys;

		while (x <= xe) {

			SetPixel(hdc, x, y, color);

			if (d < 0) {

				y += incy;
				d += d1;
			}
			else {

				d += d2;
			}

			x++;
		}
	}
	else {

		if (ys > ye) {

			int temp = ys;
			ys = ye;
			ye = temp;

			temp = xs;
			xs = xe;
			xe = temp;
		}

		int incx;

		if (xs > xe) {

			incx = -1;
		}
		else {

			incx = 1;
		}

		int deltax = xe - xs;
		int deltay = ye - ys;

		int d = (2 * abs(deltax)) - deltay;

		int d1 = (2 * abs(deltax)) - (2 * deltay);
		int d2 = 2 * abs(deltax);


		int x = xs;
		int y = ys;

		while (y <= ye) {

			SetPixel(hdc, x, y, color);

			if (d > 0) {

				x += incx;
				d += d1;
			}
			else {

				d += d2;
			}

			y++;
		}
	}
}

//Circle

double dist(int x1, int y1, int x2, int y2) {

	int dx = x2 - x1;
	int dy = y2 - y1;

	dx *= dx;
	dy *= dy;

	double r = sqrt(dx + dy);

	return r;
}

void Draw8Points(HDC hdc, int xc, int yc, int x, int y, COLORREF color) {

	SetPixel(hdc, xc + x, yc + y, color);
	SetPixel(hdc, xc - x, yc + y, color);
	SetPixel(hdc, xc + x, yc - y, color);
	SetPixel(hdc, xc - x, yc - y, color);
	SetPixel(hdc, xc + y, yc + x, color);
	SetPixel(hdc, xc - y, yc + x, color);
	SetPixel(hdc, xc + y, yc - x, color);
	SetPixel(hdc, xc - y, yc - x, color);
}

//Cartesian

void CartesianCircle(HDC hdc, int xc, int yc, double R, COLORREF color) {

	int x = 0;
	int y = R;

	while (x <= y) {

		Draw8Points(hdc, xc, yc, round(x), round(y), color);

		x++;

		y = sqrt((R*R) - (x*x));
	}
}

//Polar

void PolarCircle(HDC hdc, int xc, int yc, double R, COLORREF color) {

	double dtheta = 1.0 / R;

	double x, y;

	for (double theta = 0; theta< M_PI / 4.0; theta += dtheta) {

		x = R*cos(theta);
		y = R*sin(theta);

		Draw8Points(hdc, xc, yc, round(x), round(y), color);
	}
}

//IterativePolar
void IterativePolarCircle(HDC hdc, int xc, int yc, double R, COLORREF color) {

	double dtheta = 1.0 / R;

	double ct = cos(dtheta);
	double st = sin(dtheta);

	double x = R;
	double y = 0;
	double x1;

	while (x >= y) {

		Draw8Points(hdc, xc, yc, round(x), round(y), color);

		x1 = (x*ct) - (y*st);
		y = (x*st) + (y*ct);

		x = x1;
	}
}

//Midpoint

void CircleMidPoint(HDC hdc, int xc, int yc, double R, COLORREF color) {

	int x = 0;
	int y = R;

	int d = 1 - R;
	int d1 = 3;
	int d2 = 5 - 2 * R;

	while (x <= y) {

		Draw8Points(hdc, xc, yc, x, y, color);

		if (d < 0) {

			d += d1;
			d2 += 2;
		}
		else {

			d += d2;
			y--;
			d2 += 4;
		}

		d1 += 2;
		x++;
	}
}

//Curves

//1st Degree

void CurveFirstDegree(HDC hdc, int xs, int ys, int xe, int ye, COLORREF color)
{
	SetPixel(hdc, xs, ys, color);
	for (double t = 0; t <= 1; t += .001)
	{
		double x = (xs*(1 - t)) + (xe*t);
		double y = (ys*(1 - t)) + (ye*t);
		SetPixel(hdc, round(x), round(y), color);
	}
}

//2nd Degree

void CurveSecondDegree(HDC hdc, int x1, int y1, int xt, int yt, int x2, int y2, COLORREF color) {

	int d = max(abs(x1 - x2), abs(y1 - y2));

	double dt = 1.0 / d;

	double alpha1 = x2 - x1 - xt;
	double alpha2 = y2 - y1 - yt;
	double beta1 = xt;
	double beta2 = yt;
	double gama1 = x1;
	double gama2 = y1;

	double x = x1;
	double y = y1;

	double t = 0;

	while (t <= 1) {

		SetPixel(hdc, round(x), round(y), color);

		x = alpha1*t*t + beta1*t + gama1;
		y = alpha2*t*t + beta2*t + gama2;

		t += dt;
	}
}

//3rd Degree

//Hermite

void Hermite(HDC hdc, int x1, int y1, int x2, int y2, int x3, int y3, int x4, int y4, COLORREF color){

	int a1, b1, c1, d1;
	a1 = 2 * x1 + x2 - 2 * x3 + x4;
	b1 = -3 * x1 - 2 * x2 + 3 * x3 - x4;
	c1 = x2;
	d1 = x1;
	int a2, b2, c2, d2;
	a2 = 2 * y1 + y2 - 2 * y3 + y4;
	b2 = -3 * y1 - 2 * y2 + 3 * y3 - y4;
	c2 = y2;
	d2 = y1;
	double x, y;
	for (double t = 0; t <= 1; t += 0.0001){

		x = a1*(t*t*t) + b1*(t*t) + c1*(t)+d1;
		y = a2*(t*t*t) + b2*(t*t) + c2*(t)+d2;
		SetPixel(hdc, round(x), round(y), color);
	}
}

//Bezier

void Bezier(HDC hdc, int x1, int y1, int x2, int y2, int x3, int y3, int x4, int y4, COLORREF color)
{

	double x, y;
	for (double t = 0; t <= 1; t += 0.0001)
	{
		x = ((1 - t)*(1 - t)*(1 - t)*x1) + (3 * (1 - t)*(1 - t)*t*x2) + (3 * (1 - t)*t*t*x3) + (t*t*t*x4);
		y = ((1 - t)*(1 - t)*(1 - t)*y1) + (3 * (1 - t)*(1 - t)*t*y2) + (3 * (1 - t)*t*t*y3) + (t*t*t*y4);
		SetPixel(hdc, round(x), round(y), color);
	}
}

//nth Degree

void Spline(HDC hdc, POINT p[], int flag , COLORREF color) {

	if (flag == 2) {

		CurveSecondDegree(hdc, p[0].x, p[0].y,2*( p[2].x- p[1].x)/2,2*( p[2].y- p[1].y)/2, p[1].x, p[1].y ,color);
	}
	else if (flag > 2) {

		Hermite(hdc, p[0].x, p[0].y, 2 * (p[1].x - p[0].x) / 2, 2 * (p[1].y - p[0].y) / 2, p[1].x, p[1].y, 2 * (p[2].x - p[0].x) / 2, 2 * (p[2].y - p[0].y) / 2, color);
	}

}

//Filling

//Convex Filling

class Edge {

public:
	int xleft;
	int xright;

	Edge() {}

	Edge(int a, int b) {

		xleft = a;
		xright = b;
	}
};

void InitEdgeTable(Edge table[], int n) {

	for (int i = 0; i < n; i++) {

		table[i].xleft = 800;
		table[i].xright = 0;
	}
}

void FromEdgeToTable(POINT p1, POINT p2, Edge * table) {

	if (p1.y == p2.y) { return; }

	if (p1.y > p2.y) {

		POINT temp = p1;
		p1 = p2;
		p2 = temp;
	}

	int y = p1.y;
	double x = p1.x;

	double mi = (p2.x - p1.x) / (double)(p2.y - p1.y);

	while (y < p2.y) {

		if (x < table[y].xleft) {

			table[y].xleft = ceil(x);
		}

		if (x > table[y].xright) {

			table[y].xright = floor(x);
		}

		x += mi;
		y++;
	}
}

void FromPolygonToTable(POINT p[], int n, Edge* table) {

	POINT v1 = p[n - 1];

	for (int i = 0; i < n; i++) {

		POINT v2 = p[i];

		FromEdgeToTable(v1, v2, table);

		v1 = p[i];
	}
}

void ConvexFilling(HDC hdc, POINT p[] , int n, COLORREF color) {

	Edge table[800];

	InitEdgeTable(table, 800);

	FromPolygonToTable(p, n, table);

	for (int i = 0; i < 800; i++) {

		if (table[i].xleft <= table[i].xright) {

			LineMidPoint(hdc, table[i].xleft, i, table[i].xright, i, color);
		}
	}
}

//Clipping 

//Point Clipping

void PointClipping(HDC hdc, int x, int y, int xleft, int xright, int ytop, int ybottom, COLORREF color) {

	if (x >= xleft && x <= xright && y >= ytop && y <= ybottom) {

		SetPixel(hdc, x, y, color);
	}

}

//Line Clipping

union OutCode {

	unsigned All : 4;
	struct { unsigned left : 1, right : 1, bottom : 1, top : 1; };
};

OutCode GetOutCode(double x, double y, int xleft, int xright, int ybottom, int ytop) {

	OutCode result;

	result.All = 0;

	if (x < xleft) {

		result.left = 1;
	}
	else if (x > xright) {

		result.right = 1;
	}

	if (y > ybottom) {

		result.bottom = 1;
	}
	else if (y < ytop) {

		result.top = 1;
	}

	return result;
}

struct point {

	double x;
	double y;
};

point Vintersect(double xs, double ys, double xe, double ye, int xedge) {

	point v;

	v.x = xedge;
	v.y = ((xedge - xs)*(ye - ys) / (double)(xe - xs)) + ys;

	return v;
}

point Hintersect(double xs, double ys, double xe, double ye, int yedge) {

	point v;

	v.y = yedge;
	v.x = ((yedge - ys)*(xe - xs) / (double)(ye - ys)) + xs;

	return v;
}

void LineClipping(HDC hdc, int xs, int ys, int xe, int ye, int xleft, int xright, int ytop, int ybottom, COLORREF color) {

	OutCode out1 = GetOutCode(xs, ys, xleft, xright, ybottom, ytop);

	OutCode out2 = GetOutCode(xe, ye, xleft, xright, ybottom, ytop);

	while (((out1.All != 0 || out2.All != 0) && ((out1.All&out2.All) == 0)))
	{
		if (out1.All) {

			if (out1.left) {

				point v = Vintersect(xs, ys, xe, ye, xleft);
				xs = v.x;
				ys = v.y;
			}
			else if (out1.right) {

				point v = Vintersect(xs, ys, xe, ye, xright);
				xs = v.x;
				ys = v.y;
			}
			else if (out1.bottom) {

				point v = Hintersect(xs, ys, xe, ye, ybottom);
				xs = v.x;
				ys = v.y;
			}
			else if (out1.top) {

				point v = Hintersect(xs, ys, xe, ye, ytop);
				xs = v.x;
				ys = v.y;
			}

			out1 = GetOutCode(xs, ys, xleft, xright, ybottom, ytop);
		}
		else if (out2.All) {

			if (out2.left) {

				point v = Vintersect(xs, ys, xe, ye, xleft);
				xe = v.x;
				ye = v.y;
			}
			else if (out2.right) {

				point v = Vintersect(xs, ys, xe, ye, xright);
				xe = v.x;
				ye = v.y;
			}
			else if (out2.bottom) {

				point v = Hintersect(xs, ys, xe, ye, ybottom);
				xe = v.x;
				ye = v.y;
			}
			else if (out2.top) {

				point v = Hintersect(xs, ys, xe, ye, ytop);
				xe = v.x;
				ye = v.y;
			}

			out2 = GetOutCode(xe, ye, xleft, xright, ybottom, ytop);
		}
	}
	if (out1.All == 0 && out2.All == 0) {

		DDA(hdc, xs, ys, xe, ye, color);
	}
}

//Circle Clipping

void CirclePointClipping(HDC hdc,int xc,int yc,double R,int x,int y,COLORREF color)
{
	double d=((x - xc)*(x - xc))+ ((y - yc)*(y - yc))-(R*R);
	if (d <= 0)
	{
		SetPixel(hdc,x,y, color);
	}
}

void CircleLineClipping(HDC hdc, int xc, int yc, double R,int xs, int ys, int xe, int ye, COLORREF color)
{
	double d1 = ((xs - xc)*(xs - xc)) + ((ys - yc)*(ys - yc)) - (R*R);
	double d2 = ((xe - xc)*(xe - xc)) + ((ye - yc)*(ye - yc)) - (R*R);
	if (d1 <= 0&&d2<=0)
	{
		DDA(hdc,xs,ys,xe,ye,color);
		return;
	}
	double m = (ye - ys) / (double)(xe - xs);
	double c = ys - m*xs;
	double A, B, C;
	A = (1 + (m*m));
	B = 2 * ((m*c) - (m*yc) - (xc));
	C = (xc*xc) + ((c - yc)*(c - yc)) - (R*R);
	if ((B*B) - (4.0 * A*C) >= 0)
	{
		double x1 = ((-1 * B) - (sqrt((B*B) - (4.0 * A*C)))) / (2.0*A);
		double y1 = (m*x1) + c;
		double x2 = ((-1 * B) + (sqrt((B*B) - (4.0 * A*C)))) / (2.0*A);
		double y2 = (m*x2) + c;
		if (d1 <= 0 && d2 > 0)
		{	if(x2>=min(xs,xe)&&x2<= max(xs, xe))
				DDA(hdc, xs, ys, x2, y2, color);
			else if (x1 >= min(xs, xe) && x1 <= max(xs, xe))
				DDA(hdc, xs, ys, x1, y1, color);
		}
		else if (d1 > 0 && d2 <= 0)
		{ 
			if (x2 >= min(xs, xe) && x2 <= max(xs, xe))
				DDA(hdc,x2, y2,xe,ye, color);
			else if (x1 >= min(xs, xe) && x1 <= max(xs, xe))
				DDA(hdc, x1, y1, xe, ye, color);
		}
		else
		{	
			if (x1 >= min(xs, xe) && x1 <= max(xs, xe))
				DDA(hdc, x1, y1, x2, y2, color);
		}
	}

}

//Saving & Loading

struct Line {

	string type;
	int xs, ys, xe, ye;
};

struct Circle {

	string type;
	int xc, yc;
	double R;
};

struct Curve {

	string type;
	POINT p[100];
	int n;
};

struct Filling {

	POINT p[5];
};

vector<Line> LineList;
vector<Circle> CircleList;
vector<Curve> CurveList;
vector<Filling> FillingList;

void clear(HDC hdc, COLORREF color)
{
	for (int i = 0;i <= 685;i++)
	{
		for (int j = 0;j <= 1365;j++)
		{
			SetPixel(hdc,j,i,color);
		}
	}
}
void BackGround(HDC hdc, COLORREF pcolor, COLORREF color)
{
	for (int i = 0;i <= 685;i++)
	{
		for (int j = 0;j <= 1365;j++)
		{
			COLORREF c=GetPixel(hdc,j,i);
			if(c==pcolor)
				SetPixel(hdc, j, i, color);
		}
	}
}

void save(COLORREF color)
{
	freopen("lines.txt", "w", stdout);
	cout<<color<<endl;
	cout<< LineList.size()<<endl;
	for (int i = 0;i < LineList.size();i++)
	{
		cout << LineList[i].type << endl;
		cout<< LineList[i].xs<<" "<<LineList[i].ys<<" "<< LineList[i].xe<<" "<< LineList[i].ye<<endl;
	}
	freopen("circles.txt", "w", stdout);
	cout << CircleList.size() << endl;
	for (int i = 0;i < CircleList.size();i++)
	{
		cout << CircleList[i].type << endl;
		cout << CircleList[i].xc << " " << CircleList[i].yc << " " << CircleList[i].R << endl;
	}
	freopen("curves.txt", "w", stdout);
	cout << CurveList.size() << endl;
	for (int i = 0;i < CurveList.size();i++)
	{
		cout << CurveList[i].type << endl;
		cout<< CurveList[i].n<<endl;
		for (int j = 0;j < CurveList[i].n;j++)
		{
			cout << CurveList[i].p[j].x << " " << CurveList[i].p[j].y << endl;
		}
		
	}
	freopen("filling.txt", "w", stdout);
	cout << FillingList.size() << endl;
	for (int i = 0;i < FillingList.size();i++)
	{
		for (int j = 0;j < 5;j++)
		{
			cout << FillingList[i].p[j].x << " " << FillingList[i].p[j].y << endl;
		}
	}
}

void load(HDC hdc)
{
	
	freopen("lines.txt", "r", stdin);
	COLORREF color;
	cin>>color;
	clear(hdc,color);
	int n;
	cin>>n;
	for (int i = 0;i < n;i++)
	{
		Line l;
		cin>>l.type;
		cin>>l.xs>>l.ys>>l.xe>>l.ye;
		LineList.push_back(l);
		if (l.type == "DDA")
		{
			DDA(hdc,l.xs,l.ys,l.xe,l.ye, RGB(1, 255, 0));
		}
		else if (l.type == "Parametric")
		{
			ParametricLine(hdc, l.xs, l.ys, l.xe, l.ye, RGB(255, 1, 0));
		}
		else if (l.type == "MidPoint")
		{
			LineMidPoint(hdc, l.xs, l.ys, l.xe, l.ye, RGB(0, 150, 200));
		}
	}
	freopen("circles.txt", "r", stdin);
	cin >> n;
	for (int i = 0;i < n;i++)
	{
		Circle c;
		cin >> c.type;
		cin >> c.xc >> c.yc >> c.R;
		CircleList.push_back(c);
		if (c.type == "Cartesian")
		{
			CartesianCircle(hdc, c.xc, c.yc, c.R, RGB(0, 255, 255));
		}
		else if (c.type == "Polar")
		{
			PolarCircle(hdc, c.xc, c.yc, c.R, RGB(0, 255, 255));
		}
		else if (c.type == "IterativePolar")
		{
			IterativePolarCircle(hdc, c.xc, c.yc, c.R, RGB(0, 255, 255));
		}
		else if (c.type == "MidPoint")
		{
			CircleMidPoint(hdc, c.xc, c.yc, c.R, RGB(0, 255, 255));
		}
	}
	freopen("curves.txt", "r", stdin);
	cin >> n;
	for (int i = 0;i < n;i++)
	{
		Curve c;
		cin >> c.type;
		cin>>c.n;
		for (int j = 0;j < c.n;j++)
		{
			cin>>c.p[j].x>>c.p[j].y;
		}
		CurveList.push_back(c);
		if (c.type == "FirstDegree")
		{
			CurveFirstDegree(hdc,c.p[0].x, c.p[0].y, c.p[1].x,c.p[1].y, RGB(0, 255, 255));
		}
		else if (c.type == "SecondDegree")
		{
			CurveSecondDegree(hdc, c.p[0].x, c.p[0].y, c.p[1].x, c.p[1].y, c.p[2].x, c.p[2].y, RGB(0, 255, 255));
		}
		else if (c.type == "Hermite")
		{
			Hermite(hdc, c.p[0].x, c.p[0].y, c.p[1].x- c.p[0].x, c.p[1].y- c.p[0].y, c.p[2].x, c.p[2].y, c.p[3].x- c.p[2].x, c.p[3].y- c.p[2].y, RGB(0, 255, 255));
		}
		else if (c.type == "Bezier")
		{
			Bezier(hdc, c.p[0].x, c.p[0].y, c.p[1].x, c.p[1].y, c.p[2].x, c.p[2].y ,c.p[3].x, c.p[3].y, RGB(0, 255, 255));
		}
		else if (c.type == "Spline")
		{
			for (int j = 2;j < c.n;j++)
			{
				POINT p[3];
				p[0]=c.p[j-2];
				p[1]=c.p[j-1];
				p[2]=c.p[j];
				Spline(hdc,p,j, RGB(0, 255, 255));
			}
			
		}
	}
	freopen("filling.txt", "r", stdin);
	cin >> n;
	for (int i = 0;i < n;i++)
	{
		Filling f;
		for (int j = 0;j < 5;j++)
		{
			cin >> f.p[j].x >> f.p[j].y;
		}
		FillingList.push_back(f);
		ConvexFilling(hdc,f.p,5, RGB(0, 255, 255));
	}
}

//
//  FUNCTION: WndProc(HWND, UINT, WPARAM, LPARAM)
//
//  PURPOSE:  Processes messages for the main window.
//
//  WM_COMMAND  - process the application menu
//  WM_PAINT    - Paint the main window
//  WM_DESTROY  - post a quit message and return
//
//
LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
	static int choice=1;
	static int flag = 0;
	static int xs, ys, xe, ye,xleft , xright, ytop , ybottom,x1,y1,x2,y2;
	static COLORREF bc=RGB(255,255,255);
	static POINT p[5];
	HDC dc;
    switch (message)
    {
    case WM_COMMAND:
        {
            int wmId = LOWORD(wParam);
            // Parse the menu selections:
            switch (wmId)
            {
			case ID_FILE_SAVE:
				save(bc);
			break;
			case ID_FILE_LOAD:
				dc=GetDC(hWnd);
				load(dc);
				ReleaseDC(hWnd, dc);
				break;
			case ID_LINE_PARAM:
				flag = 0;
				choice = 0;
				break;
			case ID_LINE_DDA:
				choice = 1;
				flag = 0;
				break;
			case ID_LINE_MIDPOINT:
				choice = 2;
				flag = 0;
				break;
			case ID_CIRCLE_PARAMETRIC:
				choice = 3;
				flag = 0;
				break;
			case ID_CIRCLE_CARTESIAN:
				choice = 4;
				flag = 0;
				break;
			case ID_CIRCLE_POLAR:
				choice = 5;
				flag = 0;
				break;
			case ID_CIRCLE_ITERATIVEPOLAR:
				choice = 6;
				flag = 0;
				break;
			case ID_CIRCLE_MIDPOINT:
				choice = 7;
				flag = 0;
				break;
			case ID_CURVES_FIRSTDEGREE:
				choice = 8;
				flag = 0;
				break;
			case ID_CURVES_SECONDDEGREE:
				choice = 9;
				flag = 0;
				break;
			case ID_THIRDDEGREE_HERMITE:
				choice = 10;
				flag = 0;
				break;
			case ID_THIRDDEGREE_BEZIER:
				choice = 11;
				flag = 0;
				break;
			case ID_CURVES_SPLINE:
				choice = 12;
				flag = 0;
				break;
			case ID_FILLING_CONVEXFILLING:
				choice = 13;
				flag = 0;
				break;
			case ID_RECTANGLECLIPPING_POINTCLIPPING:
				flag = 0;
				choice = 14;
				dc = GetDC(hWnd);
				clear(dc, bc);
				ReleaseDC(hWnd, dc);
				break;
			case ID_RECTANGLECLIPPING_LINECLIPPING:
				flag = 0;
				choice = 15;
				dc = GetDC(hWnd);
				clear(dc, bc);
				ReleaseDC(hWnd, dc);
				break;
			case ID_CIRCLECLIPPING_POINTCLIPPING:
				flag = 0;
				choice = 16;
				dc = GetDC(hWnd);
				clear(dc, bc);
				ReleaseDC(hWnd, dc);
				break;
			case ID_CIRCLECLIPPING_LINECLIPPING:
				flag = 0;
				choice = 17;
				dc = GetDC(hWnd);
				clear(dc, bc);
				ReleaseDC(hWnd, dc);
				break;
			case ID_BACKGROUND_RED:
				dc = GetDC(hWnd);
				BackGround(dc, bc, RGB(255, 0, 0));
				ReleaseDC(hWnd, dc);
				bc = RGB(255, 0, 0);
				break;
			case ID_BACKGROUND_GREEN:
				dc = GetDC(hWnd);
				BackGround(dc, bc, RGB(0, 255, 0));
				ReleaseDC(hWnd, dc);
				bc = RGB(0, 255, 0);
				break;
			case ID_BACKGROUND_BLUE:
				dc = GetDC(hWnd);
				BackGround(dc, bc, RGB(0, 0, 255));
				ReleaseDC(hWnd, dc);
				bc = RGB(0, 0, 255);
				break;
			case ID_BACKGROUND_YELLOW:
				dc = GetDC(hWnd);
				BackGround(dc, bc, RGB(255, 255, 0));
				ReleaseDC(hWnd, dc);
				bc = RGB(255, 255, 0);
				break;
			case ID_BACKGROUND_WHITE:
				dc = GetDC(hWnd);
				BackGround(dc, bc, RGB(255, 255, 255));
				ReleaseDC(hWnd, dc);
				bc = RGB(255, 255, 255);
				break;
            case IDM_ABOUT:
                DialogBox(hInst, MAKEINTRESOURCE(IDD_ABOUTBOX), hWnd, About);
                break;
			/*case ID_BACKGROUND_COLOR:
				CColorDialog dlg;
				if (dlg.DoModal() == IDOK)
				{
					COLORREF color = dlg.GetColor();

				}
				break;*/
            case IDM_EXIT:
                DestroyWindow(hWnd);
                break;
            default:
                return DefWindowProc(hWnd, message, wParam, lParam);
            }
        }
        break;

	case WM_LBUTTONDOWN:
		switch (choice)
		{
		case 0:
		{
			if (flag == 0)
			{
				flag++;
				xs = LOWORD(lParam);
				ys = HIWORD(lParam);
			}
			else
			{
				flag = 0;
				xe = LOWORD(lParam);
				ye = HIWORD(lParam);
				HDC hdc = GetDC(hWnd);
				ParametricLine(hdc, xs, ys, xe, ye, RGB(255, 1, 0));
				ReleaseDC(hWnd, hdc);

				Line l;
				l.type = "Parametric";
				l.xs = xs;
				l.ys = ys;
				l.xe = xe;
				l.ye = ye;
				LineList.push_back(l);
			}
			break;
		}
		case 1:
		{
			if (flag == 0)
			{
				flag++;
				xs = LOWORD(lParam);
				ys = HIWORD(lParam);
			}
			else
			{
				flag = 0;
				xe = LOWORD(lParam);
				ye = HIWORD(lParam);
				HDC hdc = GetDC(hWnd);
				DDA(hdc, xs, ys, xe, ye, RGB(1, 255, 0));
				ReleaseDC(hWnd, hdc);

				Line l;
				l.type = "DDA";
				l.xs = xs;
				l.ys = ys;
				l.xe = xe;
				l.ye = ye;
				LineList.push_back(l);
			}
			break;
		}
		case 2:
		{
			if (flag == 0)
			{
				flag++;
				xs = LOWORD(lParam);
				ys = HIWORD(lParam);
			}
			else
			{
				flag = 0;
				xe = LOWORD(lParam);
				ye = HIWORD(lParam);
				HDC hdc = GetDC(hWnd);
				LineMidPoint(hdc, xs, ys, xe, ye, RGB(0, 150, 200));
				ReleaseDC(hWnd, hdc);

				Line l;
				l.type = "MidPoint";
				l.xs = xs;
				l.ys = ys;
				l.xe = xe;
				l.ye = ye;
				LineList.push_back(l);
			}
			break;
		}
		case 3:
		{
			if (flag == 0)
			{
				flag++;
				xs = LOWORD(lParam);
				ys = HIWORD(lParam);
			}
			else
			{
				flag = 0;
				xe = LOWORD(lParam);
				ye = HIWORD(lParam);

				double R = dist(xs, ys, xe, ye);

				HDC hdc = GetDC(hWnd);
				//parametric circle 
				ReleaseDC(hWnd, hdc);
			}
			break;
		}
		case 4:
		{
			if (flag == 0)
			{
				flag++;
				xs = LOWORD(lParam);
				ys = HIWORD(lParam);
			}
			else
			{
				flag = 0;
				xe = LOWORD(lParam);
				ye = HIWORD(lParam);

				double R = dist(xs, ys, xe, ye);

				HDC hdc = GetDC(hWnd);
				CartesianCircle(hdc, xs, ys, R, RGB(0, 255, 255));
				ReleaseDC(hWnd, hdc);

				Circle c;
				c.type = "Cartesian";
				c.xc = xs;
				c.yc = ys;
				c.R = R;
				CircleList.push_back(c);
			}
			break;
		}
		case 5:
		{
			if (flag == 0)
			{
				flag++;
				xs = LOWORD(lParam);
				ys = HIWORD(lParam);
			}
			else
			{
				flag = 0;
				xe = LOWORD(lParam);
				ye = HIWORD(lParam);

				double R = dist(xs, ys, xe, ye);

				HDC hdc = GetDC(hWnd);
				PolarCircle(hdc, xs, ys, R, RGB(0, 255, 255));
				ReleaseDC(hWnd, hdc);

				Circle c;
				c.type = "Polar";
				c.xc = xs;
				c.yc = ys;
				c.R = R;
				CircleList.push_back(c);
			}
			break;
		}
		case 6:
		{
			if (flag == 0)
			{
				flag++;
				xs = LOWORD(lParam);
				ys = HIWORD(lParam);
			}
			else
			{
				flag = 0;
				xe = LOWORD(lParam);
				ye = HIWORD(lParam);

				double R = dist(xs, ys, xe, ye);

				HDC hdc = GetDC(hWnd);
				IterativePolarCircle(hdc, xs, ys, R, RGB(0, 255, 255));
				ReleaseDC(hWnd, hdc);

				Circle c;
				c.type = "IterativePolar";
				c.xc = xs;
				c.yc = ys;
				c.R = R;
				CircleList.push_back(c);
			}
			break;
		}
		case 7:
		{
			if (flag == 0)
			{
				flag++;
				xs = LOWORD(lParam);
				ys = HIWORD(lParam);
			}
			else
			{
				flag = 0;
				xe = LOWORD(lParam);
				ye = HIWORD(lParam);

				double R = dist(xs, ys, xe, ye);

				HDC hdc = GetDC(hWnd);
				CircleMidPoint(hdc, xs, ys, R, RGB(0, 255, 255));
				ReleaseDC(hWnd, hdc);

				Circle c;
				c.type = "MidPoint";
				c.xc = xs;
				c.yc = ys;
				c.R = R;
				CircleList.push_back(c);
			}
			break;
		}
		case 8:
		{
			if (flag == 0)
			{
				flag++;
				xs = LOWORD(lParam);
				ys = HIWORD(lParam);
			}
			else
			{
				flag = 0;
				xe = LOWORD(lParam);
				ye = HIWORD(lParam);
				HDC hdc = GetDC(hWnd);
				CurveFirstDegree(hdc, xs, ys, xe , ye , RGB(0, 255, 255));
				ReleaseDC(hWnd, hdc);

				Curve c;
				c.type = "FirstDegree";
				c.n = 2;
				c.p[0].x = xs;
				c.p[0].y = ys;
				c.p[1].x = xe;
				c.p[1].y = ye;
				CurveList.push_back(c);
			}
			break;
		}
		case 9:
		{
			if (flag == 0)
			{
				
				p[0].x = LOWORD(lParam);
				p[0].y = HIWORD(lParam);
				flag++;
			}
			else if (flag == 1) {

				p[1].x = LOWORD(lParam);
				p[1].y = HIWORD(lParam);
				flag++;
			}
			else{

				flag = 0;
				p[2].x = LOWORD(lParam);
				p[2].y = HIWORD(lParam);
				HDC hdc = GetDC(hWnd);
				CurveSecondDegree(hdc, p[0].x , p[0].y , p[1].x-p[0].x , p[1].y-p[0].y , p[2].x, p[2].y, RGB(0, 255, 255));
				ReleaseDC(hWnd, hdc);
				
				Curve c;
				c.type = "SecondDegree";
				c.n = 3;
				c.p[0].x = p[0].x;
				c.p[0].y = p[0].y;
				c.p[1].x = p[1].x - p[0].x;
				c.p[1].y = p[1].y - p[0].y;
				c.p[2].x = p[2].x;
				c.p[2].y = p[2].y;
				CurveList.push_back(c);
			}
			break;
		}
		case 10:
		{
			if (flag == 0)
			{

				p[0].x = LOWORD(lParam);
				p[0].y = HIWORD(lParam);
				flag++;
			}
			else if (flag == 1) {

				p[1].x = LOWORD(lParam);
				p[1].y = HIWORD(lParam);
				flag++;
			}
			else if (flag == 2) {

				p[2].x = LOWORD(lParam);
				p[2].y = HIWORD(lParam);
				flag++;

			}else{

				flag = 0;
				p[3].x = LOWORD(lParam);
				p[3].y = HIWORD(lParam);
				HDC hdc = GetDC(hWnd);
				Hermite(hdc, p[0].x, p[0].y, (p[1].x - p[0].x), (p[1].y - p[0].y), p[2].x, p[2].y, (p[2].x - p[3].x), (p[2].y - p[3].y), RGB(0, 255, 255));
				ReleaseDC(hWnd, hdc);

				Curve c;
				c.type = "Hermite";
				c.n = 4;
				c.p[0].x = p[0].x;
				c.p[0].y = p[0].y;
				c.p[1].x = p[1].x - p[0].x;
				c.p[1].y = p[1].y - p[0].y;
				c.p[2].x = p[2].x;
				c.p[2].y = p[2].y;
				c.p[3].x = p[2].x - p[3].x;
				c.p[3].y = p[2].y - p[3].y;
				
				CurveList.push_back(c);

			}
			break;
		}
		case 11:
		{
			if (flag == 0)
			{

				p[0].x = LOWORD(lParam);
				p[0].y = HIWORD(lParam);
				flag++;
			}
			else if (flag == 1) {

				p[1].x = LOWORD(lParam);
				p[1].y = HIWORD(lParam);
				flag++;
			}
			else if (flag == 2) {

				p[2].x = LOWORD(lParam);
				p[2].y = HIWORD(lParam);
				flag++;

			}
			else {

				flag = 0;
				p[3].x = LOWORD(lParam);
				p[3].y = HIWORD(lParam);
				HDC hdc = GetDC(hWnd);
				Bezier(hdc, p[0].x, p[0].y, p[1].x , p[1].y, p[2].x, p[2].y, p[3].x,  p[3].y, RGB(0, 255, 255));
				ReleaseDC(hWnd, hdc);

				Curve c;
				c.type = "Bezier";
				c.n = 4;
				c.p[0].x = p[0].x;
				c.p[0].y = p[0].y;
				c.p[1].x = p[1].x;
				c.p[1].y = p[1].y;
				c.p[2].x = p[2].x;
				c.p[2].y = p[2].y;
				c.p[3].x = p[3].x;
				c.p[3].y = p[3].y;

				CurveList.push_back(c);
			}
			break;
		}
		case 12:
		{
			if (flag == 0)
			{

				p[0].x = LOWORD(lParam);
				p[0].y = HIWORD(lParam);
				flag++;
			}
			else if (flag == 1) {

				p[1].x = LOWORD(lParam);
				p[1].y = HIWORD(lParam);
				flag++;

				Curve c;
				c.type = "Spline";
				c.n = 2;
				c.p[0].x = p[0].x;
				c.p[0].y = p[0].y;
				c.p[1].x = p[1].x;
				c.p[1].y = p[1].y;

				CurveList.push_back(c);
			}
			else if (flag == 2) {

				p[2].x = LOWORD(lParam);
				p[2].y = HIWORD(lParam);

				HDC hdc = GetDC(hWnd);
				Spline(hdc, p, flag, RGB(0, 255, 255));
				ReleaseDC(hWnd, hdc);

				CurveList[CurveList.size() - 1].p[flag].x = p[2].x;
				CurveList[CurveList.size() - 1].p[flag].y = p[2].y;
				CurveList[CurveList.size() - 1].n++;

				flag++;
			}
			else {

				p[0] = p[1];
				p[1] = p[2];
				p[2].x = LOWORD(lParam);
				p[2].y = HIWORD(lParam);
				HDC hdc = GetDC(hWnd);
				Spline(hdc, p, flag, RGB(0, 255, 255));
				ReleaseDC(hWnd, hdc);

				CurveList[CurveList.size() - 1].p[flag].x = p[2].x;
				CurveList[CurveList.size() - 1].p[flag].y = p[2].y;
				CurveList[CurveList.size() - 1].n++;

				flag++;
			}
			break;
		}
		case 13:
		{
			if (flag < 4){

				p[flag].x = LOWORD(lParam);
				p[flag].y = HIWORD(lParam);
				flag++;
			}
			else{

				p[flag].x = LOWORD(lParam);
				p[flag].y = HIWORD(lParam);
				flag = 0;

				HDC hdc = GetDC(hWnd);
				ConvexFilling(hdc, p, 5, RGB(50, 255, 255));
				ReleaseDC(hWnd, hdc);

				Filling f;
				f.p[0].x = p[0].x ;
				f.p[0].y = p[0].y ;
				f.p[1].x = p[1].x ;
				f.p[1].y = p[1].y ;
				f.p[2].x = p[2].x ;
				f.p[2].y = p[2].y ;
				f.p[3].x = p[3].x ;
				f.p[3].y = p[3].y ;
				f.p[4].x = p[4].x ;
				f.p[4].y = p[4].y ;

				FillingList.push_back(f);
			}
			break;
		}
		case 14:
		{
			if (flag == 0) {

				xleft = LOWORD(lParam);
				flag++;
			}
			else if (flag == 1) {

				xright = LOWORD(lParam);
				flag++;
			}
			else if (flag == 2) {

				ytop = HIWORD(lParam);
				flag++;
			}
			else if (flag == 3) {

				ybottom = HIWORD(lParam);
				HDC hdc = GetDC(hWnd);
				DDA(hdc, xleft, ytop, xright, ytop, RGB(0, 150, 200));
				DDA(hdc, xright, ytop, xright, ybottom, RGB(0, 150, 200));
				DDA(hdc, xleft, ytop, xleft, ybottom, RGB(0, 150, 200));
				DDA(hdc, xleft, ybottom, xright, ybottom, RGB(0, 150, 200));
				ReleaseDC(hWnd, hdc);
				flag++;
			}
			else {

				xs = LOWORD(lParam);
				ys = HIWORD(lParam);
				HDC hdc = GetDC(hWnd);
				PointClipping(hdc, xs, ys, xleft, xright, ytop, ybottom, RGB(255, 0, 0));
				ReleaseDC(hWnd, hdc);
				//flag++;
			}
			break;
		}
		case 15:
		{
			if (flag == 0) {

				xleft = LOWORD(lParam);
				flag++;
			}
			else if (flag == 1) {

				xright = LOWORD(lParam);
				flag++;
			}
			else if (flag == 2) {

				ytop = HIWORD(lParam);
				flag++;
			}
			else if (flag == 3) {

				ybottom = HIWORD(lParam);
				HDC hdc = GetDC(hWnd);
				DDA(hdc, xleft, ytop, xright, ytop, RGB(0, 150, 200));
				DDA(hdc, xright, ytop, xright, ybottom, RGB(0, 150, 200));
				DDA(hdc, xleft, ytop, xleft, ybottom, RGB(0, 150, 200));
				DDA(hdc, xleft, ybottom, xright, ybottom, RGB(0, 150, 200));
				ReleaseDC(hWnd, hdc);
				flag++;
			}
			else if( flag==4 ){

				xs = LOWORD(lParam);
				ys = HIWORD(lParam);
				flag ++ ;
			}
			else {

				xe = LOWORD(lParam);
				ye = HIWORD(lParam);
				HDC hdc = GetDC(hWnd);

				LineClipping(hdc, xs, ys, xe, ye, xleft, xright, ytop, ybottom, RGB(0, 255, 255));

				ReleaseDC(hWnd, hdc);

				flag = 4;
			}
			break;
		}
		case 16:
		{
			if (flag == 0) {

				x1 = LOWORD(lParam);
				y1 = HIWORD(lParam);
				flag++;
			}
			else if (flag == 1) {

				x2 = LOWORD(lParam);
				y2 = HIWORD(lParam);
				double R = dist(x1, y1, x2, y2);
				HDC hdc = GetDC(hWnd);
				CircleMidPoint(hdc,x1,y1,R,RGB(0,255,255));
				ReleaseDC(hWnd, hdc);
				flag++;
			}
			else {

				xs = LOWORD(lParam);
				ys = HIWORD(lParam);
				double R = dist(x1, y1, x2, y2);
				HDC hdc = GetDC(hWnd);
				CirclePointClipping(hdc, x1, y1, R, xs, ys, RGB(255, 0, 0));
				ReleaseDC(hWnd, hdc);
			}
			break;
		}
		case 17:
		{
			if (flag == 0) {

				x1 = LOWORD(lParam);
				y1 = HIWORD(lParam);
				flag++;
			}
			else if (flag == 1) {

				x2 = LOWORD(lParam);
				y2 = HIWORD(lParam);
				double R = dist(x1, y1, x2, y2);
				HDC hdc = GetDC(hWnd);
				CircleMidPoint(hdc, x1, y1, R, RGB(0, 255, 255));
				ReleaseDC(hWnd, hdc);
				flag++;
			}
			else if (flag == 2) {

				xs = LOWORD(lParam);
				ys = HIWORD(lParam);
				flag++;
			}
			else {

				xe = LOWORD(lParam);
				ye = HIWORD(lParam);
				HDC hdc = GetDC(hWnd);
				double R = dist(x1, y1, x2, y2);
				CircleLineClipping(hdc, x1, y1, R, xs, ys, xe, ye,RGB(255, 0, 0));

				ReleaseDC(hWnd, hdc);

				flag = 2;
			}
			break;
		}
		}
		break;

    case WM_PAINT:
        {
            PAINTSTRUCT ps;
            HDC hdc = BeginPaint(hWnd, &ps);
            // TODO: Add any drawing code that uses hdc here...
            EndPaint(hWnd, &ps);
        }
        break;
    case WM_DESTROY:
        PostQuitMessage(0);
        break;
    default:
        return DefWindowProc(hWnd, message, wParam, lParam);
    }
    return 0;
}

// Message handler for about box.
INT_PTR CALLBACK About(HWND hDlg, UINT message, WPARAM wParam, LPARAM lParam)
{
    UNREFERENCED_PARAMETER(lParam);
    switch (message)
    {
    case WM_INITDIALOG:
        return (INT_PTR)TRUE;

    case WM_COMMAND:
        if (LOWORD(wParam) == IDOK || LOWORD(wParam) == IDCANCEL)
        {
            EndDialog(hDlg, LOWORD(wParam));
            return (INT_PTR)TRUE;
        }
        break;
    }
    return (INT_PTR)FALSE;
}
