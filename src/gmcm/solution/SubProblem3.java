package gmcm.solution;

import java.awt.Color;

import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdDraw;
import edu.princeton.cs.algs4.StdOut;

public class SubProblem3
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

	public static void print(int termNum) // 二维作图
	{
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
			StdDraw.setPenColor(Color.black);

			StdDraw.circle(bs[i].x, bs[i].y, dist[termNum][i]);
		}
		StdDraw.setPenRadius(0.01);
		StdDraw.setPenColor(Color.blue);
		StdDraw.point(xaver, yaver);
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

		for (double x = bs[m].x - r; x <= bs[m].x + r; x += delta)
			for (double y = bs[m].y - r; y <= bs[m].y + r; y += delta)
			{
				// for (double z = -r; z <= +r; z += delta)
				// {
				Terminal t = new Terminal(x, y); // 待定点
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
				// }
			}
		// Z = height(termNum, X, Y, a);
		// StdDraw.setPenRadius(0.01); // 轨迹作图
		// StdDraw.point(X, Y);
		StdOut.printf("%.2f\t%.2f\n", X, Y);
	}

	@SuppressWarnings("deprecation")
	public static void main(String[] args)
	{
		String file = "case025_input.txt";
		double[] input = In.readDoubles(file);
		init(input);
		// StdDraw.setCanvasSize();
		// StdDraw.setXscale(-1024,1024);
		// StdDraw.setYscale(-1024,1024);
		for (int i = 0; i < M; i++)
			SmartSearch(i);
		// StdOut.println(file + " done!");
	}
}