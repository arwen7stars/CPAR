import java.util.concurrent.ThreadLocalRandom;

public class Matrix {
	private int rows;
	private int cols;
	private int[][] matrix;

	public Matrix(int rows, int cols) {
		this.rows = rows;
		this.cols = cols;
		this.matrix = new int[rows][cols];
	}

	public void fillMatrix() {
		for(int i = 0; i < rows; i++)
			for(int j = 0; j < cols; j++){
				matrix[i][j] = ThreadLocalRandom.current().nextInt(1, 100 + 1);
			}
	}

	public void printMatrix() {
		for(int i = 0; i < rows; i++)
		{
			for(int j = 0; j < cols; j++)
				System.out.print(matrix[i][j] + "\t");
			System.out.print("\n");
		}
	}

	public void addVal(int row, int col, int val) {
		matrix[row][col] += val;
	}

	public int getVal(int row, int col) {
		return matrix[row][col];	
	}

	public Matrix multiplyNxN(Matrix input, int size) {
		Matrix result = new Matrix(size, size);

		for(int i = 0; i < size; i++) {
			for(int j = 0; j < size; j++) {
				for(int k = 0; k < size; k++) {
					result.addVal(i, j, matrix[i][k]*input.getVal(k, j));
				}
			}
		}

		return result;
	}

	public Matrix multiplyNxNAlternative(Matrix input, int size) {
		Matrix result = new Matrix(size, size);

		for(int i = 0; i < size; i++) {
			for(int k = 0; k < size; k++) {
				for(int j = 0; j < size; j++) {
					result.addVal(i, j, matrix[i][k]*input.getVal(k, j));
				}
			}
		}
		return result;
	}
}
