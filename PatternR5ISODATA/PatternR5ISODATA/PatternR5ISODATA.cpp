// PatternR5ISODATA.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "math.h"
#include "stdlib.h"

struct point
{
	double x, y;
};

struct ljz{
	double mm;
	int i, j;
};

struct llc{
	double x;//最大标准差
	int i;//第几分量
};

double distance(point a, point b){
	double dis;
	dis = sqrt((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y));
	return dis;
};

int K, stN, L, I, C;
point *Point1;
llc *st1max;
double stS, stC;//参数
int dd;//迭代次数
int *NN;
double *DJ, D;
int *aa;

void replace(int k, int n, int l, int i, double s, double c){
	K = k;
	stN = n;
	L = l;
	I = i;
	stS = s;
	stC = c;
};

int min(double *c)//求最近距离
{
	double  oo = 100000;
	int uu = 0;
	for (int i = 0; i<C; i++)
	{
		if (c[i]<oo){
			oo = c[i];
			uu = i;
		}
	}
	return uu;
}

void Update(point *a, point *c, int N)//A为样本点，B为聚类中心，C用来修改聚类中心,N为样本点总个数
{
	if (aa)
	{
		delete[]aa;
	}
	aa = new int[(C + 10) * 50];
	D = 0;
	if (NN)
	{
		delete[]NN;
	}
	NN = new int[C + 10];//每类的数量
	if (DJ)
	{
		delete[]DJ;
	}
	DJ = new double[C + 10];//每一类中样本到聚类中心的距离平均值
	if (st1max)
	{
		delete[]st1max;
	}
	st1max = new llc[C + 10];
	int zxz;

	point *st1 = new point[C];//每一类别中样本与聚类中心距离的标准差向量 
	int j;
	for (int i = 0; i<C; i++)
	{
		NN[i] = 0;
		c[i].x = 0;
		c[i].y = 0;
		DJ[i] = 0;
		st1[i].x = 0;
		st1[i].y = 0;
	}

	double *dis = new double[C];//为样本到聚类中心的距离

	for (int i = 0; i<C; i++)
	{
		for (int j = 0; j<50; j++)
		{
			aa[C*i * 50 + j] = 0;
		}
	}
	for (int i = 0; i<N; i++)
	{
		for (j = 0; j<C; j++)
		{
			dis[j] = distance(Point1[j], a[i]);
		}
		zxz = min(dis);//为距离I样本最近的聚类中心点号
		NN[zxz]++;
		aa[zxz * 50 + NN[zxz] - 1] = i;
		c[zxz].x += a[i].x;
		c[zxz].y += a[i].y;
	}
	for (int i = 0; i<C; i++)
	{
		Point1[i].x = c[i].x / NN[i];
		Point1[i].y = c[i].y / NN[i];
	}
	dd++;//迭代次数加1
	for (int i = 0; i<C; i++)
	{
		for (j = 0; j<NN[i]; j++)
		{
			DJ[i] += distance(Point1[i], a[aa[i * 50 + j]]);
			st1[i].x += (a[aa[i * 50 + j]].x - Point1[i].x)*(a[aa[i * 50 + j]].x - Point1[i].x);//计算每一类别中样本与聚类中心距离标准差
			st1[i].y += (a[aa[i * 50 + j]].y - Point1[i].y)*(a[aa[i * 50 + j]].y - Point1[i].y);
		}
		st1[i].x = sqrt(st1[i].x / NN[i]);
		st1[i].y = sqrt(st1[i].y / NN[i]);
		if (st1[i].x >= st1[i].y)//如果条件满足，则记录最大的标准差st1[i].x,并标记为第一分量
		{
			st1max[i].x = st1[i].x;
			st1max[i].i = 1;
		}
		else{
			st1max[i].x = st1[i].y;
			st1max[i].i = 2;
		}
		D += DJ[i];
		DJ[i] = DJ[i] / NN[i];
	}
	D = D / N;
	delete[]dis;
	delete[]st1;
}

void merge(point *a)
{
	ljz *m_a = new ljz[C*(C - 1) / 2];
	for (int yuy = 0; yuy<C*(C - 1) / 2; yuy++)
	{
		m_a[yuy].mm = 0;
	}
	int UUU;
	UUU = 0;
	int i, j;
	for (i = 0; i<C - 1; i++)//计算类间距离对于类间距离满足小于规定值得两个类，记录类的序号，并记录类间距离
	{
		for (j = i + 1; j<C; j++)
		{
			if (distance(Point1[i], Point1[j])<stC)
			{
				UUU++;
				m_a[UUU - 1].mm = distance(Point1[i], Point1[j]);
				m_a[UUU - 1].i = i;
				m_a[UUU - 1].j = j;
			}
		}
	}
	ljz bbbb;
	for (i = 0; i<UUU - 1; i++)//给记录的值排序
	{
		for (j = i + 1; j<UUU; j++)
		{
			if (m_a[i].mm>m_a[j].mm)
			{
				bbbb.mm = m_a[i].mm;
				bbbb.i = m_a[i].i;
				bbbb.j = m_a[i].j;
				m_a[i].mm = m_a[j].mm;
				m_a[i].i = m_a[j].i;
				m_a[i].j = m_a[j].j;
				m_a[j].mm = bbbb.mm;
				m_a[j].i = bbbb.i;
				m_a[j].j = bbbb.j;
			}
		}
	}
	int oooo;
	point YYY, *m_r;

	for (i = 0; i<UUU; i++)//从小到大合并类
	{
		if (m_a[i].mm != 0)
		{
			YYY.x = 0;
			YYY.y = 0;
			oooo = NN[m_a[i].i];
			NN[m_a[i].i] += NN[m_a[i].j];
			for (j = 0; j<NN[m_a[i].j]; j++){
				aa[m_a[i].i * 50 + oooo + j] = aa[m_a[i].j * 50 + j];
			}
			for (int TT = 0; TT<NN[m_a[i].i]; TT++)
			{
				YYY.x += a[aa[m_a[i].i * 50 + TT]].x;
				YYY.y += a[aa[m_a[i].i * 50 + TT]].y;
			}
			Point1[m_a[i].i].x = YYY.x / NN[m_a[i].i];
			Point1[m_a[i].i].y = YYY.y / NN[m_a[i].i];
			C--;
			m_r = Point1;
			Point1 = new point[C];
			oooo = 0;
			for (int TT = 0; TT<C + 1; TT++)
			{
				if (TT != m_a[i].j)
				{
					Point1[oooo].x = m_r[TT].x;
					Point1[oooo].y = m_r[TT].y;
					oooo++;
				}
			}
			delete[]m_r;
			int *NN1;
			NN1 = NN;
			NN = new int[C];
			oooo = 0;
			for (int TT = 0; TT<C + 1; TT++)
			{
				if (TT != m_a[i].j)
				{
					NN[oooo] = NN1[TT];
					oooo++;
				}
			}
			delete[]NN1;

			int *aa1;
			aa1 = aa;
			aa = new int[C * 50];
			oooo = 0;
			for (int TT = 0; TT<C + 1; TT++)
			{
				if (TT != m_a[i].j)
				{
					for (int UIO = 0; UIO<50; UIO++)
					{
						aa[oooo * 50 + UIO] = aa1[TT * 50 + UIO];
					}
					oooo++;
				}
			}
			delete[]aa1;
		}
		for (int uyt = i + 1; uyt<UUU; uyt++)
		{
			if (m_a[uyt].i == m_a[i].i || m_a[uyt].j == m_a[i].j)
			{
				m_a[uyt].mm = 0;
			}
		}
	}
	delete[]m_a;
}

int division()//分裂采用M=0.5
{
	int i, j, pp = 0, OOOO = C;
	point *point1, *po;
	for (i = 0; i<OOOO; i++)
	{
		if (st1max[i].x>stS)
		{
			if ((DJ[i]>D&&NN[i]>2 * (stN + 1)) || C <= (K / 2))
			{
				C++;
				pp++;
				point1 = new point[C];
				if (st1max[i].i == 1)
				{
					for (j = 0; j<C; j++)
					{
						if (j == C - 1)
						{
							point1[j].x = Point1[i].x + 0.5*st1max[i].x;
							point1[j].y = Point1[i].y;
						}
						else
						{
							if (j == i)
							{
								point1[j].x = Point1[i].x - 0.5*st1max[i].x;
								point1[j].y = Point1[i].y;
							}
							else
							{
								point1[j].x = Point1[i].x;
								point1[j].y = Point1[i].y;
							}
						}
					}
				}
				else
				{
					for (j = 0; j<C; j++)
					{
						if (j == C - 1)
						{
							point1[j].x = Point1[i].x;
							point1[j].y = Point1[i].y + 0.5*st1max[i].x;
						}
						else
						{
							if (j == i)
							{
								point1[j].x = Point1[i].x;
								point1[j].y = Point1[i].y - 0.5*st1max[i].x;
							}
							else
							{
								point1[j].x = Point1[j].x;
								point1[j].y = Point1[j].y;
							}
						}
					}
				}
				po = Point1;
				Point1 = point1;
				point1 = po;

				delete[]point1;
			}
		}
	}
	if (pp == 0)
	{
		return 0;
	}
	else
	{
		return 1;
	}
}


void main( )
{
	int N;//样本总数
	int ouou = 0;
	dd = 0;//迭代次数初始值为0
	point *Point0;
	printf("请输入样本总数:\n");
	scanf("%d", &N);
	Point0 = new point[N];
	for (int i = 0; i<N; i++)
	{
		printf("请输入第%d个点的坐标:\n", i + 1);
		scanf("%lf,%lf", &Point0[i].x, &Point0[i].y);
	}
	printf("请输入初始聚类中心个数:\n");
	scanf("%d", &C);
	int *a = new int[C];
	printf("请输入初始聚类中心序号:\n");
	int j;
	for (j = 0; j<C; j++)
	{
		scanf("%d", &a[j]);
	}
	Point1 = new point[C];
	for (j = 0; j<C; j++)
	{
		Point1[j] = Point0[a[j] - 1];
	}
	delete[]a;
	point *Point2;
	int j1, h, b, f;
	double c, v;
	printf("请输入初始参数格式为(K,stN,stS,stC,L,I)\n");
	scanf("%d,%d,%lf,%lf,%d,%d", &j1, &h, &c, &v, &b, &f);
	replace(j1, h, b, f, c, v);
	int ooo;
	
	while (1)
	{
		ooo = 0;
		Point2 = new point[C];

		Update(Point0, Point2, N);
		if (dd == I)
		{
			merge(Point0);
			for (int i = 0; i<C; i++)
			{
				printf("属于第%d类的点有:\n", i + 1);
				for (j = 0; j<NN[i]; j++)
				{
					if (j == NN[i] - 1)
					{
						printf("(%lf,%lf)\n", Point0[aa[i * 50 + j]].x, Point0[aa[i * 50 + j]].y);
					}
					else
					{
						printf("(%lf,%lf)\t", Point0[aa[i * 50 + j]].x, Point0[aa[i * 50 + j]].y);
					}
				}
			}
			break;
		}
		else
		{
			if (C <= K / 2)
			{
				ooo = division();
				if (ooo == 0)
				{
					merge(Point0);
				}
			}
			else
			{
				ouou = dd % 2;
				if (ouou == 0 || C >= 2 * K)
				{
					merge(Point0);
				}
				else
				{
					ooo = division();
					if (ooo == 0)
					{
						merge(Point0);
					}
				}
			}
		}
		delete[]Point2;
	}
	delete[]Point1;
	delete[]Point0;
	delete[]aa;
	delete[]NN;
	delete[]DJ;
	delete[]st1max;
	system("pause");
}