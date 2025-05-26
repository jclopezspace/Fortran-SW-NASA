#!/bin/bash

# On Linux you need to run this first to ensure the script has execute permissions:
# chmod +x _compileAndRunDotnet.sh

# Compiles the .NET project
dotnet build

# Run the compiled .NET project
dotnet run --project ModernFortran.csproj