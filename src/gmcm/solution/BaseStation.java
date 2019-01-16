package gmcm.solution;

public class BaseStation
{
	double x;
	double y;
	double z; // 坐标
	int num; // 编号
	int dimention; // 维数
	int include; // 叠合度

	boolean isNLOS;

	public BaseStation(double a, double b, double c)
	{
		x = a;
		y = b;
		z = c;
		dimention = 3;
		isNLOS = false;
	}

	public BaseStation(double a, double b)
	{
		x = a;
		y = b;
		z = 0.0;
		dimention = 2;
		isNLOS = false;
	}
}
