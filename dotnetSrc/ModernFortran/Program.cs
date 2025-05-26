
Console.WriteLine("_________________________________________________________");
Console.WriteLine();
Console.WriteLine("This program was written and compiled by Steve Winward");
Console.WriteLine("James Madison University, Mathematics Department");
Console.WriteLine("_________________________________________________________");
Console.WriteLine();

# region ACSI ART :)

Console.WriteLine("       _    ");
Console.WriteLine("      /_\\  ");
Console.WriteLine("     /___\\ ");
Console.WriteLine("    |=   =| ");
Console.WriteLine("    | NASA| ");
Console.WriteLine("    |     | ");
Console.WriteLine("    |     | ");
Console.WriteLine("    |     | ");
Console.WriteLine("    |     | ");
Console.WriteLine("    |     | ");
Console.WriteLine("    |     | ");
Console.WriteLine("    |     | ");
Console.WriteLine("    |     | ");
Console.WriteLine("    |     | ");
Console.WriteLine("    |     | ");
Console.WriteLine("    |     | ");
Console.WriteLine("    |     | ");
Console.WriteLine("   /|##!##|\\");
Console.WriteLine("  / |##!##| \\");
Console.WriteLine(" /  |##!##|  \\");
Console.WriteLine("|  / ^ | ^ \\  |");
Console.WriteLine("| /  ( | )  \\ |");
Console.WriteLine("|/   ( | )   \\|");
Console.WriteLine("    ((   ))");
Console.WriteLine("   ((  :  ))");
Console.WriteLine("   ((  :  ))");
Console.WriteLine("    ((   ))");
Console.WriteLine("     (( ))");
Console.WriteLine("      ( )");
Console.WriteLine("       .");
Console.WriteLine("       .");
Console.WriteLine("       .");
Console.WriteLine();

# endregion

string filename = "FILTER_IN.DAT";

// Read in the filter data
int n;
double[] s = FortranModule.ReadInFile(filename, out n);

Console.WriteLine($"The file {filename} was successfully loaded with {n} data values");
Console.WriteLine();

// Allocate memory for tridiagonal system vectors
int nMinus2 = n - 2;
double[] lVector = new double[nMinus2];
double[] uVector = new double[nMinus2];
double[] dVector = new double[nMinus2];
double[] bVector = new double[nMinus2];
double[] xInitial = new double[nMinus2];

// Initialize all vectors to 0.0
for (int i = 0; i < nMinus2; i++)
{
    lVector[i] = 0.0;
    uVector[i] = 0.0;
    dVector[i] = 0.0;
    bVector[i] = 0.0;
    xInitial[i] = 0.0;
}

// Prompt user for gamma
Console.Write("Please input a value for gamma (-1, 1): ");
string? gammaInput = Console.ReadLine();
double gamma;
while (!double.TryParse(gammaInput, out gamma))
{
    Console.Write("Invalid input. Please enter a valid number for gamma (-1, 1): ");
    gammaInput = Console.ReadLine();
}

Console.WriteLine("");

// Obtain the right hand side of the equation
FortranModule.GetRHS(gamma, n, s, bVector);

// Set up the tridiagonal matrix
FortranModule.InitializeTridiagMatrix(uVector, dVector, lVector, gamma, nMinus2);

// Allocate for iterations
int[] iterations = new int[99];

Console.WriteLine("(W, Iterations)");

# region Scatter Plot of W vs the number of iterations

for (int i = 1; i <= 99; i++)
{
    double w = 1.0 + i * 0.01;
    iterations[i - 1] = FortranModule.SuccOverRelaxation(lVector, dVector, uVector, bVector, xInitial, w, nMinus2);
    for (int j = 0; j < xInitial.Length; j++) xInitial[j] = 0.0;

    Console.WriteLine($"{w:F2}, {iterations[i - 1]}");
}

int maxIterations = 0;
for (int i = 0; i < iterations.Length; i++)
{
    if (iterations[i] > maxIterations)
        maxIterations = iterations[i];
}

int plotHeight = 20; // Number of rows for the y-axis
int plotWidth = iterations.Length; // One column per w value

// Scale function for y-axis
int ScaleY(int value)
{
    return plotHeight - 1 - (int)Math.Round((double)value / maxIterations * (plotHeight - 1));
}

// Build plot
char[,] plot = new char[plotHeight, plotWidth];
for (int y = 0; y < plotHeight; y++)
    for (int x = 0; x < plotWidth; x++)
        plot[y, x] = ' ';

// Plot points
for (int x = 0; x < plotWidth; x++)
{
    int y = ScaleY(iterations[x]);
    plot[y, x] = '*';
}

// Print y-axis labels and plot
Console.WriteLine();
for (int y = 0; y < plotHeight; y++)
{
    int iterLabel = (int)Math.Round(maxIterations * (double)(plotHeight - 1 - y) / (plotHeight - 1));
    Console.Write($"{iterLabel,4} | ");
    for (int x = 0; x < plotWidth; x++)
        Console.Write(plot[y, x]);
    Console.WriteLine();
}

// Print x-axis
Console.Write("     +");
for (int x = 0; x < plotWidth; x++)
    Console.Write('-');
Console.WriteLine();

// Print x-axis labels (w values at start, middle, end)
double wStart = 1.01;
double wEnd = 2.0;
int mid = plotWidth / 2;
Console.Write("      ");
for (int x = 0; x < plotWidth; x++)
{
    if (x == 0)
        Console.Write($"{wStart:F2}");
    else if (x == mid)
        Console.Write($"{wStart + (wEnd - wStart) / 2:F2}".PadLeft(mid - 2));
    else if (x == plotWidth - 1)
        Console.Write($"{wEnd:F2}".PadLeft(plotWidth - mid - 4));
    else
        Console.Write(" ");
}
Console.WriteLine();

# endregion