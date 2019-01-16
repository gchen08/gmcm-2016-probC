package gmcm.solution;

import java.util.Arrays;

import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdOut;
import edu.princeton.cs.algs4.StdRandom;

public class SubProblem2
{

	static final double speed = 3.0 * Math.pow(10, 8); // 传播速度

	static int N; // 基站个数
	static int M; // 终端个数
	static int dimen; // 维数

	static BaseStation[] bs; // 基站组
	static double[][] TOA; // TOA矩阵
	static double[][] dist; // 距离矩阵
	static double[][] bsDist; // 基站间距离

	static final int threshold = 4; // 最少基站数目
	static boolean[] valid; // 使用基站标识
	static double[] bsLen; // 基站模长

	public static void init(double[] data) // 初始化
	{
		N = (int) data[0]; // 读入基站数目
		M = (int) data[1]; // 读入终端数目
		dimen = (int) data[2]; // 读入维数
		bs = new BaseStation[N];
		TOA = new double[M][N];
		dist = new double[M][N];
		bsDist = new double[N][N];
		bsLen = new double[N];
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
		for (int i = 0; i < N; i++)
		{
			bsDist[i][i] = Double.MAX_VALUE;
			for (int j = i + 1; j < N; j++)
			{
				bsDist[i][j] = distance(bs[i], bs[j]);
				bsDist[j][i] = bsDist[i][j]; // 对称性
			}
		}
		double bsX = 0.0;
		double bsY = 0.0;
		double bsZ = 0.0;
		for (int i = 0; i < N; i++)
		{
			bsX += bs[i].x;
			bsY += bs[i].y;
			bsZ += bs[i].z;
		}
		bsX /= N;
		bsY /= N;
		bsZ /= N;
		BaseStation bsCenter = new BaseStation(bsX, bsY, bsZ);
		for (int i = 0; i < N; i++)
		{
			bsLen[i] = distance(bs[i], bsCenter);
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
		int m = -1;
		for (int i = 0; i < N; i++)
			if (valid[i])
			{
				m = i;
				break;
			}
		for (int i = m + 1; i < N; i++)
		{
			if ((valid[i]) && (dist[termNum][i] < dist[termNum][m]))
				m = i;
		}
		return m;
		
//		int m = 0;
//		for (int i = 1; i < N; i++)
//		{
//			if (dist[termNum][i] < dist[termNum][m])
//				m = i;
//		}
//		return m;
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
			if (valid[i])
				r[i] = Math.sqrt((x - bs[i].x) * (x - bs[i].x) + (y - bs[i].y) * (y - bs[i].y));
		double[] li = new double[N];
		double minH = Double.MAX_VALUE;
		for (int i = 0; i < N; i++)
			if (valid[i])
				minH = Math.min(minH, bs[i].z);
		// minH /= threshold;
		// StdOut.println(minH);
		int round = 1;
		double averZ = 0.0;
		for (int j = 0; j < round; j++)
		{
			for (int i = 0; i < N; i++)
			{
				if (valid[i])
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

	public static void Combine(int current)
	{
		int diff = current - threshold;
		if (diff > 0)
		{
			double[] dist = new double[N];
			for (int i = 0; i < N; i++)
				if (valid[i])
				{
					dist[i] = bsLen[i];
				}
			Arrays.sort(dist);
			for (int i = 0; i < diff; i++)
				for (int j = 0; j < N; j++)
					if ((valid[j]) && (bsLen[j] == dist[i]))
					{
						valid[j] = false;
						current--;
					}
		}
		else if (diff < 0)
		{
			double[] dist = new double[N];
			for (int i = 0; i < N; i++)
				if (!valid[i])
				{
					dist[i] = bsLen[i];
				}
			diff *= -1;
			Arrays.sort(dist);
			for (int i = 0; i < diff; i++)
				for (int j = 0; j < N; j++)
					if ((!valid[j]) && (bsLen[j] == dist[dist.length - i - 1]))
					{
						valid[j] = true;
						current++;
					}
		}
	}

	public static void SelectSearch(int termNum) // 筛选基站定位
	{
		valid = new boolean[N]; // 初始化

		filter(termNum, dist[termNum]);
		int current = 0;
		for (int i = 0; i < N; i++) // NLOS筛选
			if (!bs[i].isNLOS)
			{
				current++;
				valid[i] = true;
			}

		if (current != threshold)
		{
			int include[] = new int[N]; // 叠合度判断
			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
					if (dist[termNum][i] > bsDist[i][j])
						include[i]++;
			if (current > threshold)
			{
				loop: for (int i = 0; i < N; i++)
					if ((valid[i]) && (include[i] > threshold))
					{
						valid[i] = false;
						current--;
						if (current == threshold)
						{
							break loop;
						}
					}
			}
			else
			{
				loop: for (int i = 0; i < N; i++)
					if ((include[i] <= threshold) && !valid[i])
					{
						valid[i] = true;
						current++;
						if (current == threshold)
						{
							break loop;
						}
					}
			}
			Combine(current);
		}
		for (int i = 0; i < N; i++)
			if (valid[i])
				StdOut.printf("%d\t", bs[i].num);
		// StdOut.println();

		int m = minBaseStation(termNum); // 搜索中心
		double r = dist[termNum][m]; // 搜索半径
		double delta = 1; // 搜索步长
		double denominator = 0.0; // 参数分母
		for (int i = 0; i < N; i++)
		{
			if (valid[i])
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
						if (valid[i])
							numerator += distance(t, bs[i]);
					}
					double lemda = numerator / denominator; // 参数
					double tmp = 0.0; // 最小化
					for (int i = 0; i < N; i++)
					{
						if (valid[i])
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
		String file = "case020_input.txt";
		double[] input = In.readDoubles(file);
		init(input);
		for (int i = 0; i < M; i++)
			SelectSearch(i);
		StdOut.println(file + " done!");
	}
}