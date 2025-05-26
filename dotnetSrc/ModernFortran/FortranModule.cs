public static class FortranModule
{
    public static double ConditionNum(double gamma)
    {
        return (1 + Math.Abs(gamma)) / (1 - Math.Abs(gamma));
    }

    public static void GetRHS(double gamma, int n, double[] sUnfiltered, double[] rhs)
    {
        double a = (1.0 + gamma) / 2.0;
        double aOver2 = a / 2.0;

        rhs[0] = a * sUnfiltered[1] + aOver2 * sUnfiltered[2] + (aOver2 - gamma / 2.0) * sUnfiltered[0];

        for (int i = 1; i < n - 3; i++)
        {
            rhs[i] = aOver2 * sUnfiltered[i] + a * sUnfiltered[i + 1] + aOver2 * sUnfiltered[i + 2];
        }
        int last = n - 3;
        rhs[last] = aOver2 * sUnfiltered[last] + a * sUnfiltered[last + 1] + (aOver2 - gamma / 2.0) * sUnfiltered[last + 2];
    }

    public static void InitializeTridiagMatrix(double[] uVector, double[] dVector, double[] lVector, double gamma, int n)
    {
        uVector[0] = gamma / 2.0;
        dVector[0] = 1.0;
        lVector[0] = 0.0;
        for (int i = 1; i < n - 1; i++)
        {
            lVector[i] = gamma / 2.0;
            uVector[i] = lVector[i];
            dVector[i] = 1.0;
        }
        lVector[n - 1] = gamma / 2.0;
        dVector[n - 1] = 1.0;
        uVector[n - 1] = 0.0;
    }

    public static double[] ReadInFile(string filename, out int n)
    {
        using (var reader = new StreamReader(filename))
        {
            string? line = reader.ReadLine();
            if (!int.TryParse(line, out n))
                throw new Exception("Could not read n from file.");
            n++;
            double[] s = new double[n];
            for (int i = 0; i < n; i++)
            {
                if (!double.TryParse(reader.ReadLine(), out s[i]))
                    throw new Exception($"Could not read value #{i} from file.");
            }
            return s;
        }
    }

    public static int SuccOverRelaxation(double[] lVector, double[] dVector, double[] uVector, double[] bVector, double[] xInitial, double w, int n)
    {
        double[] xNew = new double[n];
        double[] residual = new double[n];
        double[] tempVector = new double[n];
        int maxIterations = 1000;
        double tolerance = Math.Sqrt(0.000000000000001);
        double bVectorNorm = NormInfinity(bVector, n);
        double toleranceProduct = tolerance * bVectorNorm;
        int iterations = 0;

        for (int i = 0; i < maxIterations; i++)
        {
            iterations = i + 1;

            TriDiagMatrixTimesVector(lVector, dVector, uVector, xInitial, tempVector, n);
            for (int j = 0; j < n; j++)
                residual[j] = bVector[j] - tempVector[j];

            double residualNorm = NormInfinity(residual, n);
            if (residualNorm < toleranceProduct)
            {
                break;
            }

            // First element
            xNew[0] = w * (bVector[0] - uVector[0] * xInitial[1]);
            xNew[0] = (1.0 - w) * xInitial[0] + xNew[0];

            // Middle elements
            for (int j = 1; j < n - 1; j++)
            {
                xNew[j] = w * (bVector[j] - lVector[j] * xNew[j - 1] - uVector[j] * xInitial[j + 1]);
                xNew[j] = (1.0 - w) * xInitial[j] + xNew[j];
            }

            // Last element
            xNew[n - 1] = w * (bVector[n - 1] - lVector[n - 1] * xNew[n - 2]);
            xNew[n - 1] = (1.0 - w) * xInitial[n - 1] + xNew[n - 1];

            for (int j = 0; j < n; j++)
                xInitial[j] = xNew[j];
        }

        return iterations;
    }

    public static void TriDiagMatrixTimesVector(double[] lVector, double[] dVector, double[] uVector, double[] xVector, double[] resultantVector, int n)
    {
        resultantVector[0] = dVector[0] * xVector[0] + uVector[0] * xVector[1];
        for (int i = 1; i < n - 1; i++)
        {
            resultantVector[i] = lVector[i] * xVector[i - 1] + dVector[i] * xVector[i] + uVector[i] * xVector[i + 1];
        }
        resultantVector[n - 1] = lVector[n - 1] * xVector[n - 2] + dVector[n - 1] * xVector[n - 1];
    }

    public static double NormInfinity(double[] x, int n)
    {
        double max = Math.Abs(x[0]);
        for (int i = 1; i < n; i++)
        {
            if (Math.Abs(x[i]) > max)
                max = Math.Abs(x[i]);
        }
        return max;
    }
}