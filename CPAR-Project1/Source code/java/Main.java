import java.util.concurrent.ThreadLocalRandom;

public class Main {
	public static void main(String[] args) {
		int operation, size;	
		boolean loop = true;
	
		do {
	    		if (args.length == 2) {
				try {
					operation = Integer.parseInt(args[0]);
					size = Integer.parseInt(args[1]);				
				} catch(NumberFormatException e) {
					System.err.println("Error! Arguments must be integers.");
					break;		
				}

				multiplyMatricesNxN(size, operation);
				break;
			} else {
				System.err.println("usage: Main <operation> <size>");
				System.err.println("\toperation 1 : performance of multiplication of NxN matrices (regular algorthim)");
				System.err.println("\toperation 2 : performance of multiplication of NxN matrices (alternative algorithm)");
				System.err.println("\toperation 3 : ratio of performance between different algorithms");
				System.err.println("\tsize : matrix size");
				break;	
			}
		} while(loop);		

	}
	
	public static void multiplyMatricesNxN(int matrix_size, int operation) {
		System.out.println("MULTIPLICATION OF MATRIXES " + matrix_size + "X" + matrix_size);
		double firstElapsedTime, secondElapsedTime, ratio;
		Matrix m1 = new Matrix(matrix_size, matrix_size);
		Matrix m2 = new Matrix(matrix_size, matrix_size);

		m1.fillMatrix();
		m2.fillMatrix();

		switch(operation) {
			case 1:
				System.out.println("\n\tOp1 : performance of multiplication of NxN matrices (regular algorithm)\n");
				checkPerformance(m1, m2, matrix_size, true);
				break;
			case 2:
				System.out.println("\n\tOp2 : performance of multiplication of NxN matrices (alternative algorithm)\n");
				checkPerformance(m1, m2, matrix_size, false);
				break;
			case 3:
				System.out.println("\n\tOp3 : ratio of performance between different algorithms\n");
				firstElapsedTime = checkPerformance(m1, m2, matrix_size, true);
				secondElapsedTime = checkPerformance(m1, m2, matrix_size, false);

				ratio = (double)(firstElapsedTime)/(double)(secondElapsedTime);
				System.out.println("\nRatio between first and second cycle's execution times: " + ratio + "\n");
				break;		
			default:
				break;
		}		
	}

	public static long checkPerformance(Matrix m1, Matrix m2, int matrix_size, boolean regular) {
		long startTime, stopTime, elapsedTime;
		Matrix result;

		startTime = System.currentTimeMillis();
		if (regular)
			result = m1.multiplyNxN(m2, matrix_size);
		else result = m1.multiplyNxNAlternative(m2, matrix_size);
		stopTime = System.currentTimeMillis();
		elapsedTime = stopTime - startTime;
		
		if (regular)
			System.out.println("Regular multiplication execution time: " + elapsedTime + " miliseconds");
		else System.out.println("Alternative multiplication execution time: " + elapsedTime + " miliseconds");

		// m1.printMatrix();
		// m2.printMatrix();
		// result.printMatrix();

		return elapsedTime;
	}
}
