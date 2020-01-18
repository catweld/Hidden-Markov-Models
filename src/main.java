import java.io.File;
import java.util.Arrays;
import java.util.Scanner;

public class main {
	public static double[][] matrixmultiply(double[][] matrix1, double[][] matrix2) {
		double[][] result = new double[matrix1.length][matrix2[0].length];
		for (int i = 0; i < matrix1.length; ++i) {
			for (int j = 0; j < matrix2.length; ++j) {
				for (int k = 0; k < matrix2[0].length; ++k) {
					result[i][j] += matrix1[i][k] * matrix2[k][j];
				}
			}
		}
		return result;
	}

	public static double[] createArray(String exe) {
		String[] split = exe.split(" ");
		double[] array = new double[split.length];
		for (int i = 0; i < array.length; i++) {
			float value = Float.parseFloat(split[i]);
			array[i] = value;
		}
		return array;
	}

	public static double[][] createMat(double[] array, int rows, int cols) {
		int indexArray = 2;
		double[][] matrix = new double[rows][cols];
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				matrix[i][j] = array[indexArray];
				indexArray++;
			}
		}
		return matrix;
	}

	public static double[][] transpose(double[][] matrix) {
		double[][] tpose = new double[matrix[0].length][matrix.length];
		for (int i = 0; i < matrix[0].length; i++) {
			for (int j = 0; j < matrix.length; j++) {
				tpose[i][j] = matrix[j][i];
			}
		}
		return tpose;
	}

	public static double[][] elementmultiply(double[][] matrix1, double[][] matrix2, int col1, int col2) {
		double[][] result = new double[matrix1.length][matrix1[0].length];
		for (int i = 0; i < matrix1.length; i++) {
			result[i][0] = matrix1[i][col1] * matrix2[i][col2];
		}
		return result;
	}

	public static double[][] initialAlpha(double[][] matrixPi, double[][] matrixB, int[] arrayObsInt) {
		double[][] alpha1 = new double[matrixPi[0].length][1];
		for (int j = 0; j < matrixPi[0].length; ++j) {
			alpha1[j][0] = matrixPi[0][j] * matrixB[j][arrayObsInt[1]];
		}
		return alpha1;
	}

	public static double cfactorAlpha1(double[][] alpha) {
		double c0 = 0;
		for (int i = 0; i < alpha.length; i++) {
			c0 = c0 + alpha[i][0];
		}
		return c0;
	}

	public static double[][] scaleAlpha(double[][] alpha, double c) {
		double[][] scaledalpha = new double[alpha.length][1];
		c = 1 / c;
		for (int i = 0; i < alpha.length; i++) {
			scaledalpha[i][0] = c * alpha[i][0];
		}
		return scaledalpha;
	}

	public static double[][] alphapass(double[][] alpha1, double c0, double[][] matrixA, double[][] matrixB, int acol,
			int[] arrayObsInt) {
		double[][] alpha = new double[alpha1.length][1];
		double[][] alphastore = new double[alpha1.length + 1][arrayObsInt.length - 1];
		double[][] alphaA = new double[alpha1.length][1];

		double ct;
		for (int f = 0; f < alpha1.length; f++) {
			alphastore[f][0] = alpha1[f][0];
			alphastore[alphastore.length - 1][0] = 1 / c0;
		}
		alpha = alpha1;
		for (int t = 1; t < arrayObsInt.length - 1; t++) {
			ct = 0;
			alpha = transpose(alpha);
			alphaA = matrixmultiply(alpha, matrixA);
			alphaA = transpose(alphaA);
			alpha = elementmultiply(alphaA, matrixB, acol, arrayObsInt[t + 1]);
			for (int i = 0; i < alpha.length; i++) {
				ct = ct + alpha[i][0];
			}
			alphastore[alphastore.length - 1][t] = 1 / ct;
			alpha = scaleAlpha(alpha, ct);
			for (int j = 0; j < alpha.length; j++) {
				alphastore[j][t] = alpha[j][0];
			}
		}
		return alphastore;
	}

	public static double[][] initialBeta(double[][] matrixA, double[][] cT, int[] arrayObsInt) {
		double[][] betaT = new double[matrixA.length][1];
		for (int j = 0; j < matrixA.length; j++) {
			betaT[j][0] = cT[arrayObsInt.length - 2][0];
		}
		return betaT;
	}

	public static double[][] betapass(double[][] betaT, double[][] matrixA, double[][] matrixB, int[] arrayObsInt,
			double[][] c) {
		double[][] betastore = new double[betaT.length][arrayObsInt.length - 1];
		double[][] beta = new double[matrixA.length][arrayObsInt.length - 1];
		double ct = 0;
		for (int i = 0; i < beta.length; i++) {
			for (int j = 0; j < beta[0].length; j++) {
				beta[i][j] = 0;
			}
		}
		for (int i = 0; i < betaT.length; i++) {
			betastore[i][arrayObsInt.length - 2] = betaT[i][0];
		}
		for (int t = arrayObsInt[0] - 2; t >= 0; t--) {
			for (int i = 0; i < matrixA.length; i++) {
				beta[i][t] = 0;
				for (int j = 0; j < matrixA[0].length; j++) {
					beta[i][t] += matrixA[i][j] * matrixB[j][arrayObsInt[t + 2]] * betastore[j][t + 1];
				}
				ct = c[t][0];
				beta[i][t] = ct * beta[i][t];
				// System.out.print(beta[i][t] + " ");
				betastore[i][t] = beta[i][t];
			}
		}
		return betastore;
	}

	public static double[][] gammapass(double[][] alphastore, double[][] betastore, double[][] matrixA,
			double[][] matrixB, int[] arrayObsInt) {
		double[][] gammastore = new double[matrixA.length][arrayObsInt.length-1];
		double[][] gamma = new double[matrixA.length][matrixA.length];
		for (int t = 0; t < arrayObsInt.length - 2; t++) {
			for (int i = 0; i < matrixA.length; i++) {
				gammastore[i][t] = 0;
				for (int j = 0; j < matrixA[0].length; j++) {
					gamma[i][j] = alphastore[i][t] * matrixA[i][j] * matrixB[j][arrayObsInt[t + 2]]*betastore[j][t + 1];
					gammastore[i][t] = gammastore[i][t] + gamma[i][j];
				}
			}
		}
		for (int i = 0; i < matrixA.length; i++) {
			gammastore[i][arrayObsInt.length - 2] = alphastore[i][arrayObsInt.length - 2];
		}
		return gammastore;
	}

	public static double[][][] gammadi(double[][] alphastore, double[][] betastore, double[][] matrixA,
			double[][] matrixB, int[] arrayObsInt) {
		double[][][] gamma = new double[matrixA.length][matrixA.length][arrayObsInt.length - 1];
		for (int t = 0; t < arrayObsInt.length - 2; t++) {
			for (int i = 0; i < matrixA.length; i++) {
				for (int j = 0; j < matrixA[0].length; j++) {
					gamma[i][j][t] = alphastore[i][t] * matrixA[i][j] * matrixB[j][arrayObsInt[t + 2]]* betastore[j][t + 1];
				}
			}
		}
		return gamma;
	}

	public static double[][] reestimatepi(double[][] gamma) {
		double[][] pinew = new double[1][gamma.length];
		for (int i = 0; i < gamma.length; i++) {
			pinew[0][i] = gamma[i][0];
		}
		return pinew;
	}

	public static double[][] reestimatea(double[][] matrixA, double[][] gamma, double[][][] gammadi,
			int[] arrayObsInt) {
		double[][] anew = new double[matrixA.length][matrixA[0].length];
		for (int i = 0; i < gamma.length; i++) {
			double denom = 0;
			for (int t = 0; t < arrayObsInt.length - 1; t++) {
				denom = denom + gamma[i][t];
			}
			for (int j = 0; j < matrixA.length; j++) {
				double numer = 0;
				for (int t = 0; t < arrayObsInt.length - 1; t++) {
					numer = numer + gammadi[i][j][t];
				}
				anew[i][j] = numer / denom;
			}
		}
		return anew;
	}

	public static double[][] reestimateb(double[][] matrixB, double[][] gamma, int[] arrayObsInt) {
		double[][] bnew = new double[matrixB.length][matrixB[0].length];
		for (int i = 0; i < matrixB.length; i++) {
			double denom = 0;
			for (int t = 0; t < arrayObsInt.length - 1; t++) {
				denom = denom + gamma[i][t];
			}
			for (int j = 0; j < matrixB[0].length; j++) {
				double numer = 0;
				for (int t = 0; t < arrayObsInt.length - 1; t++) {
					if (arrayObsInt[t+1] == j) {
						numer = numer + gamma[i][t];
					}
				}
				bnew[i][j] = numer / denom;
			}
		}
		return bnew;
	}

	public static double logprob(double[][] c) {
		double logc = 0;
		for (int i = 0; i < c.length; i++) {
			logc = logc + Math.log(c[i][0]);
		}
		logc = -logc;
		return logc;
	}

	public static void main(String[] args) throws Exception {
		File input = new File("C:\\Users\\Catherine\\Desktop\\Java\\HMM3 - Estimate Model\\hmm_c_N1000.txt");
		Scanner in = new Scanner(input);
		String A = in.nextLine();
		String B = in.nextLine();
		String P = in.nextLine();
		String O = in.nextLine();
		int arows, acols, brows, bcols, prows, pcols;

		// Create arrays from strings
		double[] arrayA = createArray(A);
		double[] arrayB = createArray(B);
		double[] arrayPi = createArray(P);
		double[] arrayObs = createArray(O);

		// Assign rows and columns
		arows = (int) Math.round(arrayA[0]);
		acols = (int) Math.round(arrayA[1]);
		brows = (int) Math.round(arrayB[0]);
		bcols = (int) Math.round(arrayB[1]);
		prows = (int) Math.round(arrayPi[0]);
		pcols = (int) Math.round(arrayPi[1]);

		// Create matrices A,B,Pi
		double[][] matrixA = createMat(arrayA, arows, acols);
		double[][] matrixB = createMat(arrayB, brows, bcols);
		double[][] matrixPi = createMat(arrayPi, prows, pcols);

		// Initialize by multiplying (element-wise) transition matrix with current (initial) estimate of states
		int[] arrayObsInt = new int[arrayObs.length];
		for (int i = 0; i < arrayObsInt.length; i++) {
			arrayObsInt[i] = (int) Math.round(arrayObs[i]);
		}

		//initialize variables
		double[][] alpha1 = new double[pcols][1];
		double[][] alpha = new double[alpha1.length + 1][arrayObsInt.length - 1];
		double[][] scaledalpha = new double[alpha1.length][1];
		double[][] c = new double[alpha[0].length][1];
		double[][] betaT = new double[matrixPi[0].length][1];
		double[][] beta = new double[alpha.length][1];
		double[][] gamma = new double[matrixA.length][arrayObsInt.length - 1];
		double[][][] gammadi = new double[matrixA.length][matrixA.length][arrayObsInt.length - 1];

		int maxiters = 1000;
		double oldlogprob = -10000000;
		double logprob = 0;
		int iterations = 0;

		for (int iters = 0; (iters < maxiters) & (logprob >= oldlogprob); iters++) {
			if (iters == 0) {
				oldlogprob = -100000;
			}
			else {oldlogprob = logprob;}
			iterations = iterations+1;
			alpha1 = initialAlpha(matrixPi, matrixB, arrayObsInt);
			double c0 = cfactorAlpha1(alpha1);
			scaledalpha = scaleAlpha(alpha1, c0);

			alpha = alphapass(scaledalpha, c0, matrixA, matrixB, 0, arrayObsInt);
			
			//retrieve c array
			for (int i = 0; i < alpha[0].length; i++) {
				c[i][0] = alpha[alpha.length - 1][i];
			}

			//initialize beta
			betaT = initialBeta(matrixA, c, arrayObsInt);

			//run the beta pass
			beta = betapass(betaT, matrixA, matrixB, arrayObsInt, c);

			//run the gamma pass
			gammadi = gammadi(alpha, beta, matrixA, matrixB, arrayObsInt);
			gamma = gammapass(alpha, beta, matrixA, matrixB, arrayObsInt);

			//reestimate pi
			matrixPi = reestimatepi(gamma);

			//reestimate a
			matrixA = reestimatea(matrixA, gamma, gammadi, arrayObsInt);

			//reestimate b
			matrixB = reestimateb(matrixB, gamma, arrayObsInt);

			//compute log
			logprob = logprob(c);
			
		}
		System.out.print(iterations);
		System.out.print(matrixA.length + " " + matrixA[0].length + " ");
		for (int i = 0; i < matrixA.length; i++) {
			for (int j = 0; j < matrixA[0].length; j++) {
				System.out.print(matrixA[i][j] + " ");
			}
		}
		System.out.println();
		System.out.print(matrixB.length + " " + matrixB[0].length + " ");
		for (int i = 0; i < matrixB.length; i++) {
			for (int j = 0; j < matrixB[0].length; j++) {
				System.out.print(matrixB[i][j] + " ");
			}
		}
		

	}
}