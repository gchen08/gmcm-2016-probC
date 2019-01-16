package gmcm.solution;

import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdOut;

public class Compare
{

	static final double speed = 3.0 * Math.pow(10, 8); // 传播速度

	static int N; // 基站个数
	static int M; // 终端个数
	static int dimen; // 维数

	static BaseStation[] bs; // 基站组
	static double[][] TOA; // TOA矩阵
	static double[][] dist; // 距离矩阵

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

	public static void fit(double x, double y, double z)
	{
		double[] r = new double[N];
		for (int i = 0; i < N; i++)
			r[i] = Math.sqrt(
					(x - bs[i].x) * (x - bs[i].x) + (y - bs[i].y) * (y - bs[i].y) + (z - bs[i].z) * (z - bs[i].z));
		final double vdt = 5.0;
		// dr
		StdOut.print("[");
		for (int i = 0; i < N - 1; i++)
		{
			StdOut.print(vdt * r[i] + " ");
		}
		StdOut.println(vdt * r[N - 1] + "]");
		// dx
		StdOut.print("[");
		for (int i = 0; i < N - 1; i++)
		{
			StdOut.print(x - bs[i].x + " ");
		}
		StdOut.println(x - bs[N - 1].x + "]");
		// dy
		StdOut.print("[");
		for (int i = 0; i < N - 1; i++)
		{
			StdOut.print(y - bs[i].y + " ");
		}
		StdOut.println(y - bs[N - 1].y + "]");
		// dz
		StdOut.print("[");
		for (int i = 0; i < N - 1; i++)
		{
			StdOut.print(z - bs[i].z + " ");
		}
		StdOut.println(z - bs[N - 1].z + "]");
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

	public static void calc(int termNum, double[] ansData, double[] ourData)
	{
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

		int p = (termNum - 1) * 3;
		Terminal T = new Terminal(ourData[p], ourData[p + 1], ourData[p + 2]);
		double[] TDist = new double[N];
		for (int i = 0; i < N; i++)
			TDist[i] = distance(T, bs[i]);
		double delta = 0.0;
		for (int i = 0; i < N; i++)
			delta += 1.0 / TDist[i];
		delta *= 3.0 / N;
		StdOut.println(delta);
		StdOut.printf("%.2f\t%.2f\t%.2f\n", delta * (T.x - bsX), delta * (T.y - bsY), delta * (T.z - bsZ));
		StdOut.printf("%.2f\t%.2f\t%.2f\n", ansData[p], ansData[p + 1], ansData[p + 2]);
		StdOut.printf("%.2f\t%.2f\t%.2f\n", ourData[p], ourData[p + 1], ourData[p + 2]);
	}

	public static void diff(double[] A, double[] B)
	{
		if (A.length == B.length)
		{
			for (int i = 0; i < A.length; i += 3)
			{
				Terminal a = new Terminal(A[i], A[i + 1], A[i + 2]);
				BaseStation b = new BaseStation(B[i], B[i + 1], B[i + 2]);
				StdOut.printf("%.3f\n", distance(a, b));
			}
		}
	}

	public static void format(double[] data) // 格式化数据
	{
		for (int i = 0; i < data.length; i+=3)
		{
			StdOut.printf("%.2f\t%.2f\t%.2f\n", data[i], data[i + 1], data[i + 2]);
		}
	}

	@SuppressWarnings("deprecation")
	public static void main(String[] args)
	{
		// String inFile = "sample_case002_input.txt";
		// double[] data = In.readDoubles(inFile);
		// init(data);
//		String ansFile = "sample_case005_ans.txt";
		String ourFile = "output_case_020.txt";
//		double[] ansData = In.readDoubles(ansFile);
		double[] ourData = In.readDoubles(ourFile);
		// int termNum = 5; // 测试点
		// calc(termNum, ansData, ourData);
		// fit(x, y, z); //多元回归分析
		// diff(ansData, ourData); //误差测量
		format(ourData);
	}

}
