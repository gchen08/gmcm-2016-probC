package gmcm.solution;

public class Terminal
{
	double x;
	double y;
	double z;
	int dimention;
	boolean valid;
	
	public Terminal(double a, double b, double c)
	{
		x = a;
		y = b;
		z = c;
		dimention = 3;
	}
	
	public Terminal(double a, double b)
	{
		x = a;
		y = b;
		z = 0.0;
		dimention = 2;
	}
}
