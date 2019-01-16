package gmcm.solution;

import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdOut;
import edu.princeton.cs.algs4.StdRandom;

public class SubProblem4
{

	static final double speed = 3.0 * Math.pow(10, 8); // 传播速度

	static int N; // 基站个数
	static int M; // 终端个数
	static int dimen; // 维数

	static BaseStation[] bs; // 基站组
	static double[][] TOA; // TOA矩阵
	static double[][] dist; // 距离矩阵

	static final int threshold = 4; // 最少基站数目
	static boolean[][] valid; // 通信有效性
	static int[] connectNum; // 连接数
	static boolean[] CanbeLocated; // 被定位

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

	public static void PreTreat() // 预处理过程
	{
		// StdOut.print(M + " " + N + " ");
		valid = new boolean[M][N];
		connectNum = new int[M];
		CanbeLocated = new boolean[M];
		for (int i = 0; i < M; i++)
			for (int j = 0; j < N; j++)
				if (dist[i][j] <= 200)
				{
					valid[i][j] = true;
					connectNum[i]++;
				}
		// int countT = 0; //统计量
		// int countC = 0;
		for (int i = 0; i < M; i++)
			if (connectNum[i] >= threshold)
			{
				// countT++;
				// countC += connectNum[i];
				CanbeLocated[i] = true;
			}
		// StdOut.println(countT + " " + countC);
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
		int m = -1;
		for (int i = 0; i < N; i++)
			if (valid[termNum][i])
			{
				m = i;
				break;
			}
		for (int i = m + 1; i < N; i++)
		{
			if ((valid[termNum][i]) && (dist[termNum][i] < dist[termNum][m]))
				m = i;
		}
		return m;
	}

	public static void filter(int termNum, double[] tDist) // 初步判断基站是否NLOS
	{
		int[] bsCount = new int[N];
		for (int i = 0; i < N; i++) // 任意三个基站进行组合
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
			if (bsCount[i] > 0)
			{
				bs[i].isNLOS = true;
			}
	}

	public static double height(int termNum, double x, double y)
	{
		double[] r = new double[N]; // 平面距离
		for (int i = 0; i < N; i++)
			if (valid[termNum][i])
				r[i] = Math.sqrt((x - bs[i].x) * (x - bs[i].x) + (y - bs[i].y) * (y - bs[i].y));
		double[] li = new double[N];
		double minH = Double.MAX_VALUE;
		for (int i = 0; i < N; i++)
			if (valid[termNum][i])
				minH = Math.min(minH, bs[i].z);
		// minH /= threshold;
		// StdOut.println(minH);
		int round = 1;
		double averZ = 0.0;
		for (int j = 0; j < round; j++)
		{
			for (int i = 0; i < N; i++)
			{
				if (valid[termNum][i])
				{
					double lB = ((bs[i].z - minH) * (bs[i].z - minH) + r[i] * r[i])
							/ (dist[termNum][i] * dist[termNum][i]);
					double uB = (bs[i].z * bs[i].z + r[i] * r[i]) / (dist[termNum][i] * dist[termNum][i]);
					li[i] = StdRandom.uniform(lB, uB);
					averZ += bs[i].z - Math.sqrt(li[i] * dist[termNum][i] * dist[termNum][i] - r[i] * r[i]);
				}
			}
		}
		averZ /= threshold * round;
		return averZ;
	}

	public static void validSearch(int termNum) // 对可定位终端进行定位
	{
		if (!CanbeLocated[termNum])
		{
			StdOut.printf("%s\t%s\t%s\n", "xnan", "ynan", "znan");
		}
		else
		{
			int m = minBaseStation(termNum); // 搜索中心
			double r = dist[termNum][m]; // 搜索半径
			double delta = 1; // 搜索步长
			double denominator = 0.0; // 参数分母
//			StdOut.printf("终端：%d\t所连接基站：", termNum+1);
			for (int i = 0; i < N; i++)
				if (valid[termNum][i])
				{
//					StdOut.printf("%d\t", bs[i].num);
					denominator += dist[termNum][i];
			}
//			StdOut.println();
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
							if (valid[termNum][i])
								numerator += distance(t, bs[i]);
						}
						double lemda = numerator / denominator; // 参数
						double tmp = 0.0; // 最小化
						for (int i = 0; i < N; i++)
						{
							if (valid[termNum][i])
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
	}

	@SuppressWarnings("deprecation")
	public static void main(String[] args)
	{
		String file = "case030_input.txt";
		// StdOut.print(file + " ");
		double[] input = In.readDoubles(file);
		init(input);
		PreTreat();
		for (int i = 0; i < M; i++)
			validSearch(i);
//		StdOut.println(file + " done!");
	}
}