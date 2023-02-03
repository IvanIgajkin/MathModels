// See https://aka.ms/new-console-template for more information

using MathNet.Numerics.LinearAlgebra;

var years = new []{1902, 1907, 1912, 1917};
var population = new []{1_174_700, 1_345_700, 1_617_157, 1_854_400};

Console.WriteLine($"Population in 1910: MathSquares - {SolveByMinSquares(1910)}\tLagranzhPolinom - {SolveByLagrangePolynomial(1910)}");
Console.WriteLine($"Population in 1916: MathSquares - {SolveByMinSquares(1916)}\tLagranzhPolinom - {SolveByLagrangePolynomial(1916)}");

double SolveByMinSquares(int x) => 
	MinSquares(years.Select(x => x - 1900), population.Select(y => y / 1E6))(x - 1900);

double SolveByLagrangePolynomial(int x) => 
	LagrangePolynomial(years.Select(x => x - 1900), population.Select(y => y / 1E6))(x - 1900);

Func<int, double> MinSquares<TX, TY>(IEnumerable<TX> xData, IEnumerable<TY> yData)
{
	var x = xData.Select(xi => Convert.ToDouble(xi));
	var y = yData.Select(yi => Convert.ToDouble(yi));
	
	int n = x.Count();

	Vector<double> FindCoefficients()
	{
		var A = Matrix<double>.Build.DenseOfRows(Enumerable.Range(0, n)
			.Select(i => Enumerable.Range(0, n)
				.Select(j => x.Select(xi => Math.Pow(xi, n + 2 - (i + j))).Sum())
			));
	
		var b = Vector<double>.Build.DenseOfEnumerable(Enumerable.Range(0, n)
			.Select(degree => y.Select((yi, i) => yi * Math.Pow(x.ElementAt(i), n - 1 - degree)).Sum()));

		if (A == null || b == null)
		{
			return Vector<double>.Build.DenseOfEnumerable(Enumerable.Empty<double>());
		}
		
		return A.Solve(b);
	}

	return arg => FindCoefficients()
		.Select((c, i) => c * Math.Pow(arg, n - 1 - i))
		.Sum();
}

Func<int, double> LagrangePolynomial<TX, TY>(IEnumerable<TX> xData, IEnumerable<TY> yData)
{
	var x = xData.Select(xi => Convert.ToDouble(xi));
	var y = yData.Select(yi => Convert.ToDouble(yi));

	return arg => y.Select((yi, i) => yi * x
			.Select((xj, j) => i == j ? 1.0 : (arg - xj) / (x.ElementAt(i) - xj))
			.Aggregate(1.0, (xl, xr) => xl * xr))
		.Sum();
}
