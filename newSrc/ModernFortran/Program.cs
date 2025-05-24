
Console.WriteLine("_________________________________________________________");
Console.WriteLine();
Console.WriteLine("This program was written and compiled by Steve Winward");
Console.WriteLine("James Madison University, Mathematics Department");
Console.WriteLine("_________________________________________________________");
Console.WriteLine();

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

// Obtain the right hand side of the equation
FortranModule.GetRHS(gamma, n, s, bVector);

// Set up the tridiagonal matrix
FortranModule.InitializeTridiagMatrix(uVector, dVector, lVector, gamma, nMinus2);

// Allocate for iterations
int[] iterations = new int[99];

using (var writer = new StreamWriter("iterations1.dat"))
{
    for (int i = 1; i <= 99; i++)
    {
        double w = 1.0 + i * 0.01;
        Console.WriteLine($"W= {w:F4}");
        iterations[i - 1] = FortranModule.SuccOverRelaxation(lVector, dVector, uVector, bVector, xInitial, w, nMinus2);
        for (int j = 0; j < xInitial.Length; j++) xInitial[j] = 0.0;
        Console.WriteLine(iterations[i - 1]);
        writer.WriteLine(iterations[i - 1]);
    }
}