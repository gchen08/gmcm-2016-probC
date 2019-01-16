package gmcm.solution;

import java.awt.Color;

import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdDraw;
import edu.princeton.cs.algs4.StdOut;
import edu.princeton.cs.algs4.StdRandom;

public class SubProblem1
{

	static final double speed = 3.0 * Math.pow(10, 8); // 传播速度

	static int N; // 基站个数
	static int M; // 终端个数
	static int dimen; // 维数

	static BaseStation[] bs; // 基站组
	static double[][] TOA; // TOA矩阵
	static double[][] dist; // 距离矩阵

	static double XMIN; // 点范围
	static double XMAX;
	static double YMIN;
	static double YMAX;
	static double ZMIN;
	static double ZMAX;

	static Terminal ans = new Terminal(-21.19, 4.48, 1.48); // 准确位置

	public static void init(double[] data) // 初始化
	{
		N = (int) data[0]; // 读入基站数目
		M = (int) data[1]; // 读入终端数目
		dimen = (int) data[2]; // 读入维数
		bs = new BaseStation[N];
		TOA = new double[M][N];
		dist = new double[M][N];
		int count = 3;
		if (dimen == 3) // 三维情形
		{
			for (int i = 0; i < N; i++)
			{
				bs[i] = new BaseStation(data[count++], data[count++], data[count++]);
				bs[i].num = i + 1;
			}
		}
		else // 二维
		{
			for (int i = 0; i < N; i++)
			{
				bs[i] = new BaseStation(data[count++], data[count++]);
				bs[i].num = i + 1;
			}
		}
		for (int i = 0; i < M; i++)
			for (int j = 0; j < N; j++)
			{
				TOA[i][j] = data[count++];
				dist[i][j] = TOA[i][j] * speed;
			}
	}

	public static double distance(Terminal A, BaseStation B) // 终端A到基站B的距离
	{
		if (A.dimention == 3)
		{
			return Math.sqrt((A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y) + (A.z - B.z) * (A.z - B.z));
		}
		else
		{
			return Math.sqrt((A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y));
		}
	}

	public static double distance(BaseStation A, BaseStation B) // 基站A到基站B的距离
	{
		if (A.dimention == 3)
		{
			return Math.sqrt((A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y) + (A.z - B.z) * (A.z - B.z));
		}
		else
		{
			return Math.sqrt((A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y));
		}

	}

	public static boolean isValid(Terminal P, double[] Dist) // 判断终端P位置是否准确
	{
		double eps = 0.5; // 噪声误差
		for (int i = 0; i < N; i++)
		{
			if (distance(P, bs[i]) > Dist[i] + eps)
				return false;
		}
		return true;
	}

	public static int Com(int n) // C(n, 3)组合数
	{
		return (n) * (n - 1) * (n - 2) / 6;
	}

	public static int minBaseStation(int termNum) // 选取最近基站
	{
		int m = 0;
		for (int i = 1; i < N; i++)
		{
			if (dist[termNum][i] < dist[termNum][m])
				m = i;
		}
		return m;
	}

	public static void filter(int termNum, double[] tDist) // 判断基站是否NLOS
	{
		int[] bsCount = new int[N];
		for (int i = 0; i < N; i++)
			for (int j = i + 1; j < N; j++)
				for (int k = j + 1; k < N; k++)
				{
					double Pij = distance(bs[i], bs[j]);
					double Pik = distance(bs[i], bs[k]);
					double Pkj = distance(bs[k], bs[j]);
					if (tDist[i] > tDist[j] + Pij)
					{
						bsCount[i]++;
					}
					if (tDist[i] > tDist[k] + Pik)
					{
						bsCount[i]++;
					}
					if (tDist[j] > tDist[i] + Pij)
					{
						bsCount[j]++;
					}
					if (tDist[j] > tDist[k] + Pkj)
					{
						bsCount[j]++;
					}
					if (tDist[k] > tDist[i] + Pik)
					{
						bsCount[k]++;
					}
					if (tDist[k] > tDist[j] + Pkj)
					{
						bsCount[k]++;
					}
				}
		for (int i = 0; i < N; i++)
		{
			// StdOut.println(bsCount[i]);
			if (bsCount[i] > 0)
			{
				bs[i].isNLOS = true;
			}
			else
			{
				bs[i].isNLOS = false;
			}
		}
	}

	public static void BruteForceSearch(int termNum) // 搜索可能位置
	{
		int m = minBaseStation(termNum); // 搜索中心
		double r = dist[termNum][m]; // 搜索半径
		double delta = 0.5; // 搜索步长
		Terminal[] tl = new Terminal[100000000]; // 待筛选终端
		XMIN = Double.MAX_VALUE;
		YMIN = Double.MAX_VALUE;
		ZMIN = Double.MAX_VALUE;
		XMAX = Double.MIN_VALUE;
		YMAX = Double.MIN_VALUE;
		ZMAX = Double.MIN_VALUE;
		int count = -1; // 合法点数
		double xaver = 0.0; // 平均点
		double yaver = 0.0;
		double zaver = 0.0;
		for (double x = bs[m].x - r; x <= bs[m].x + r; x += delta) // 暴力搜索
			for (double y = bs[m].y - r; y <= bs[m].y + r; y += delta)
				for (double z = bs[m].z - r; z <= bs[m].z + r; z += delta)
				{
					if (isValid(new Terminal(x, y, z), dist[termNum]))
					{
						count++;
						tl[count] = new Terminal(x, y, z);
						tl[count].valid = true;
						xaver += x;
						yaver += y;
						zaver += z;
						if (x < XMIN)
							XMIN = x;
						if (x > XMAX)
							XMAX = x;
						if (y < YMIN)
							YMIN = y;
						if (y > YMAX)
							YMAX = y;
						if (z < ZMIN)
							ZMIN = z;
						if (z > ZMAX)
							ZMAX = z;
					}
				}
		// StdOut.println(count);
		xaver /= count;
		yaver /= count;
		zaver /= count;
		Terminal averP = new Terminal(xaver, yaver, zaver); // 平均点
		StdDraw.setCanvasSize(); // 作图
		StdDraw.setXscale(XMIN, XMAX);
		StdDraw.setYscale(YMIN, YMAX);
		StdOut.println(ans.x + " " + ans.y + " " + ans.z); // 真实点位
		StdOut.println(averP.x + " " + averP.y + " " + averP.z);
		double[] newDist = new double[N]; // 通过基站关系进行距离修正
		for (int i = 0; i < N; i++)
			newDist[i] = dist[termNum][i];
		filter(termNum, dist[termNum]);
		boolean flag = true;
		while (flag)
		{
			for (int i = 0; i < N; i++)
			{
				double lemda = 1.0;
				if (bs[i].isNLOS)
				{
					lemda = 1 - distance(averP, bs[i]) / newDist[i];
				}
				newDist[i] = (1 - lemda) * distance(averP, bs[i]) + (lemda) * newDist[i];
			}
			filter(termNum, newDist);
			flag = false;
			for (int i = 0; i < N; i++)
				if (bs[i].isNLOS)
				{
					flag = true;
					break;
				}
		}
		xaver = 0.0;
		yaver = 0.0;
		zaver = 0.0;
		double number = 0.0;
		for (int i = 0; i < count; i++)
		{
			if (isValid(tl[i], newDist))
			{
				xaver += tl[i].x;
				yaver += tl[i].y;
				zaver += tl[i].z;
				number++;
			}
		}
		xaver /= number;
		yaver /= number;
		zaver /= number;
		averP.x = xaver;
		averP.y = yaver;
		averP.z = zaver;
		// triBSLoc(termNum, newDist); //修正位置后 任三个基站定位作图
		StdOut.printf("%.2f\t%.2f\t%.2f\n", averP.x, averP.y, averP.z);
	}

	public static void triBSLoc(int termNum, double[] Dist) // 三点基站定位
	{
		for (int i = 0; i < N; i++)
			for (int j = i + 1; j < N; j++)
				for (int k = j + 1; k < N; k++)
				{
					double Ri = Dist[i]; // 到基站i的距离
					double Rj = Dist[j]; // 到基站j的距离
					double Rk = Dist[k]; // 到基站k的距离
					double ai = bs[i].x;
					double bi = bs[i].y;
					double ci = bs[i].z;
					double aj = bs[j].x;
					double bj = bs[j].y;
					double cj = bs[j].z;
					double ak = bs[k].x;
					double bk = bs[k].y;
					double ck = bs[k].z;
					double aij = ai - aj;
					double aik = ai - ak;
					double bij = bi - bj;
					double bik = bi - bk;
					double cij = ci - cj;
					double cik = ci - ck;
					double gij = 0.5
							* ((ai * ai + bi * bi + ci * ci) - (aj * aj + bj * bj + cj * cj) - (Ri * Ri - Rj * Rj));
					double gik = 0.5
							* ((ai * ai + bi * bi + ci * ci) - (ak * ak + bk * bk + ck * ck) - (Ri * Ri - Rk * Rk));
					double delta = 0.1; // 搜索步长
					double eps = 1E-8;
					loop: for (double z = ZMIN; z <= ZMAX; z += delta)
					{
						double x = (-(bij * gik - bik * gij) + (bij * cik - bik * cij) * z) / (aij * bik - aik * bij);
						double y = ((aij * gik - aik * gij) - (aij * cik - aik * cij) * z) / (aij * bik - aik * bij);
						if ((x - ai) * (x - ai) + (y - bi) * (y - bi) + (z - ci) * (z - ci) - Ri * Ri < eps)
						{
							if (isValid(new Terminal(x, y, z), Dist))
							{
								StdDraw.point(x, y);
								break loop;
							}
						}
					}
					loop: for (double z = ZMAX; z > ZMIN; z -= delta)
					{
						double x = (-(bij * gik - bik * gij) + (bij * cik - bik * cij) * z) / (aij * bik - aik * bij);
						double y = ((aij * gik - aik * gij) - (aij * cik - aik * cij) * z) / (aij * bik - aik * bij);
						if ((x - ai) * (x - ai) + (y - bi) * (y - bi) + (z - ci) * (z - ci) - Ri * Ri < eps)
						{
							if (isValid(new Terminal(x, y, z), Dist))
							{
								StdDraw.point(x, y);
								break loop;
							}
						}
					}
				}
		// StdOut.println(Count);
	}

	public static void print(int termNum) // 二维作图
	{
		filter(termNum, dist[termNum]);
		double xmin = Double.MAX_VALUE;
		double ymin = Double.MAX_VALUE;
		double xmax = Double.MIN_VALUE;
		double ymax = Double.MIN_VALUE;
		double xaver = 0.0;
		double yaver = 0.0;
		for (int i = 0; i < N; i++)
		{
			xmin = Math.min(xmin, bs[i].x - dist[termNum][i]);
			ymin = Math.min(ymin, bs[i].y - dist[termNum][i]);
			xmax = Math.max(xmax, bs[i].x + dist[termNum][i]);
			ymax = Math.max(ymax, bs[i].y + dist[termNum][i]);
			xaver += bs[i].x;
			yaver += bs[i].y;
		}
		xaver /= N;
		yaver /= N;
		StdDraw.setCanvasSize();
		StdDraw.setXscale(xmin, xmax);
		StdDraw.setYscale(ymin, ymax);
		for (int i = 0; i < N; i++)
		{
			if (!bs[i].isNLOS)
			{
				StdDraw.setPenColor(Color.black);
				// StdOut.println(bs[i].num + ": " + dist[termNum][i]);
			}
			else
			{
				StdDraw.setPenColor(Color.yellow);
			}
			StdDraw.circle(bs[i].x, bs[i].y, dist[termNum][i]);
		}
		StdDraw.setPenRadius(0.01);
		StdDraw.setPenColor(Color.red);
		StdDraw.point(ans.x, ans.y); // 实际位置
		StdDraw.setPenColor(Color.blue);
		StdDraw.point(xaver, yaver);
	}

	public static void compare(int termNum) // 准确位置偏移
	{
		for (int i = 0; i < N; i++)
		{
			double diff = distance(ans, bs[i]);
			StdOut.printf("%.2f\t%.2f\n", dist[termNum][i], diff);
		}
	}

	public static double height(int termNum, double x, double y)
	{
		double[] r = new double[N]; // 平面距离
		for (int i = 0; i < N; i++)
		{
			r[i] = Math.sqrt((x - bs[i].x) * (x - bs[i].x) + (y - bs[i].y) * (y - bs[i].y));
		}
		double[] li = new double[N];
		double minH = 0;
		for (int i = 0; i < N; i++)
			minH += bs[i].z;
		minH /= N;
		int round = 1;
		double averZ = 0.0;
		for (int j = 0; j < round; j++)
		{
			for (int i = 0; i < N; i++)
			{
				double lB = ((bs[i].z - minH) * (bs[i].z - minH) + r[i] * r[i]) / (dist[termNum][i] * dist[termNum][i]);
				double uB = (bs[i].z * bs[i].z + r[i] * r[i]) / (dist[termNum][i] * dist[termNum][i]);
				li[i] = StdRandom.uniform(lB, uB);
				averZ += bs[i].z - Math.sqrt(li[i] * dist[termNum][i] * dist[termNum][i] - r[i] * r[i]);
			}
		}
		averZ /= N * round;
		return averZ;
	}

	public static void SmartSearch(int termNum)
	{
		int m = minBaseStation(termNum); // 搜索中心
		double r = dist[termNum][m]; // 搜索半径
		double delta = 1; // 搜索步长
		double denominator = 0.0; // 参数分母
		for (int i = 0; i < N; i++)
		{
			denominator += dist[termNum][i];
		}
		double minF = Double.MAX_VALUE;
		double X = 0;
		double Y = 0;
		double Z = 0;
		for (double x = bs[m].x - r; x <= bs[m].x + r; x += delta)
			for (double y = bs[m].y - r; y <= bs[m].y + r; y += delta)
				for (double z = bs[m].z - r; z <= bs[m].z + r; z += delta)
				{
					Terminal t = new Terminal(x, y, z); // 待定点
					double numerator = 0.0; // 参数分子
					for (int i = 0; i < N; i++)
					{
						numerator += distance(t, bs[i]);
					}
					double lemda = numerator / denominator; // 参数
					double tmp = 0.0; // 最小化
					for (int i = 0; i < N; i++)
					{
						tmp += Math.pow(lemda * dist[termNum][i] - distance(t, bs[i]), 2);
					}
					if (tmp < minF)
					{
						minF = tmp;
						X = x;
						Y = y;
					}
				}
		Z = height(termNum, X, Y);
		StdOut.printf("%.2f\t%.2f\t%.2f\n", X, Y, Z);
	}

	@SuppressWarnings("deprecation")
	public static void main(String[] args)
	{
		String file = "case010_input.txt";
		double[] input = In.readDoubles(file);
		init(input);
		for (int i = 0; i < M; i++)
			SmartSearch(i);
//		StdOut.println(file + " done!");
	}
}