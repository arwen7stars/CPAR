using System;
using System.Diagnostics;
class Matrix {
	static void Main(string[] args) {
	
		if (args.Length == 1) {
			int size = Convert.ToInt32(args[0]);
			multiplyMatricesNxN(size);
		} else {
 			System.Console.WriteLine("usage: Matrix <size>");
		}
    		
  	}

	public static void multiplyMatricesNxN(int size) {
		long elapsed_time;		
		int[,] a = new int[size,size];
		int[,] b = new int[size,size];
		int[,] result = new int[size,size];

		Console.WriteLine("MULTIPLICATION OF MATRICES " + size + "X" + size);

		Random rnd = new Random();
		FillMatrix(a, rnd);
		FillMatrix(b, rnd);

		var stopwatch = new Stopwatch();
		stopwatch.Start();

		result = MultiplyMatricesNxNAux(a,b,size);

		stopwatch.Stop();
		elapsed_time = stopwatch.ElapsedMilliseconds;

		//PrintMatrix(a);
		//PrintMatrix(b);
		//PrintMatrix(result);

		Console.WriteLine("Execution time: " + elapsed_time + " miliseconds");
	}

	public static void FillMatrix(int[,] matrix, Random rnd)
	{
		for (int i = 0; i < matrix.GetLength(0); i++)
		{
			for (int j = 0; j < matrix.GetLength(1); j++)
			{
		    		matrix[i,j] = rnd.Next(1, 100);
			}
		}
	}
	
	public static void PrintMatrix(int[,] matrix)
	{
		for (int i = 0; i < matrix.GetLength(0); i++)
	    	{
			for (int j = 0; j < matrix.GetLength(1); j++)
			{
				Console.Write(matrix[i,j] + "\t");
			}
			Console.WriteLine("");
	    	}
		Console.WriteLine("");
	}

	public static int[,] MultiplyMatricesNxNAux(int[,] a, int[,] b, int size) {
		int[,] result = new int[size, size];


		for(int i = 0; i < size; i++) {
			for(int j = 0; j < size; j++) {
				for(int k = 0; k < size; k++) {
					result[i,j] += a[i,k]*b[k,j];
				}
			}
		}
		return result;
	}
}
