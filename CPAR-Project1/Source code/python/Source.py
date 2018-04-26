import random, time, sys, argparse

def multiplyMatricesNxN(size):
	a = fillMatrix(size)
	b = fillMatrix(size)
	result = initializeMatrix(size)

	print "MULTIPLICATION OF MATRIXES",size,"X",size

	start_time = time.time()
	multiplyMatricesAux(a,b,result,size)
	elapsed_time = time.time() - start_time

	#printMatrix(a, size)
	#printMatrix(b, size)
	#printMatrix(result, size)

	print "Execution time: ",elapsed_time," seconds\n"

def fillMatrix(size):
	matrix = []
	for x in range(size):
		row = []
	    	for y in range(size):
			row.append(random.randint(1,100))
	    	matrix.append(row)
	return matrix

def initializeMatrix(size):
	matrix = []
	for x in range(size):
		row = []
	    	for y in range(size):
			row.append(0)
	    	matrix.append(row)
	return matrix

def multiplyMatricesAux(a, b, result, size):
	for i in range(size):
	    	for j in range(size):
			for k in range(size):
				result[i][j] = result[i][j] + a[i][k]*b[k][j]

def printMatrix(matrix, size):
	for i in range(size):
    		for j in range(size):
        		print '{:4}'.format(matrix[i][j]),
    		print
	print
	

if __name__ == '__main__':
	if len(sys.argv) != 2:
		print "usage: ",sys.argv[0]," <size>\n"
	else:
		p = argparse.ArgumentParser()
		p.add_argument("size", type=int, help="matrix size")
		args = p.parse_args()

		multiplyMatricesNxN(args.size)



