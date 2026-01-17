#!/bin/bash
# setup_scripts/build_opencosmo.sh
# Builds the openCOSMO-RS C++ standalone binary AND Python extension
#
# Usage: ./setup_scripts/build_opencosmo.sh [path_to_cpp_repo]
# Example: ./setup_scripts/build_opencosmo.sh ~/resources/openCOSMO-RS_cpp

set -e  # Exit on error

echo "=========================================="
echo "openCOSMO-RS C++ Builder (Binary + Python)"
echo "=========================================="

# Parse arguments
if [ $# -eq 0 ]; then
    CPP_REPO="$HOME/resources/openCOSMO-RS_cpp"
    echo "No path provided, using default: $CPP_REPO"
else
    CPP_REPO="$1"
    echo "Using provided path: $CPP_REPO"
fi

# Validate path exists
if [ ! -d "$CPP_REPO" ]; then
    echo "✗ Error: Directory not found: $CPP_REPO"
    exit 1
fi

echo "✓ Repository found"

# Check for required tools
echo ""
echo "Checking build dependencies..."

if ! command -v g++ &> /dev/null; then
    echo "✗ Error: g++ not found."
    exit 1
fi
echo "✓ g++ found: $(g++ --version | head -1)"

# Check for OpenMP support
if ! echo | g++ -fopenmp -x c++ - -o /dev/null 2>/dev/null; then
    echo "⚠ Warning: OpenMP not available, building without it"
    OPENMP_FLAG=""
else
    echo "✓ OpenMP support found"
    OPENMP_FLAG="-fopenmp"
fi

# CPU instruction set
PARALLEL_FLAG="-mavx"

echo ""
echo "Build configuration:"
echo "  - Optimization: -O3"
echo "  - Parallelization: $OPENMP_FLAG"
echo "  - Instruction set: $PARALLEL_FLAG"
echo "  - Standard: C++14"

# Navigate to bindings folder
cd "$CPP_REPO/bindings"
echo ""
echo "Building in: $(pwd)"

###############################################
# 1. Build Standalone CLI Binary
###############################################

echo ""
echo "------------------------------------------"
echo "Building standalone CLI binary..."
echo "------------------------------------------"

# Remove old binary
if [ -f "openCOSMORS" ]; then
    echo "Removing old binary..."
    rm openCOSMORS
fi

# Build binary
g++ $OPENMP_FLAG $PARALLEL_FLAG -O3 -Wall -std=c++14 \
    ../code/bindings_forCLI.cpp \
    -o openCOSMORS \
    -I .. \
    -I ../eigen

# Validate binary
if [ ! -f "openCOSMORS" ]; then
    echo "✗ Binary build failed!"
    exit 1
fi

chmod +x openCOSMORS
echo "✓ CLI binary built successfully"


###############################################
# 2. Build Python Extension (.so)
###############################################

echo ""
echo "------------------------------------------"
echo "Building Python extension (.so)..."
echo "------------------------------------------"

SO_NAME="openCOSMORS$(python3-config --extension-suffix)"

# Remove old .so
if [ -f "$SO_NAME" ]; then
    echo "Removing old Python module..."
    rm "$SO_NAME"
fi

# Build .so
g++ $OPENMP_FLAG $PARALLEL_FLAG -O3 -Wall -shared -std=c++14 -fPIC \
    ../code/bindings_forPython.cpp \
    -o "$SO_NAME" \
    -I ../pybind11/include \
    -I ../eigen \
    -I ../nlohmann \
    -I /usr/include/python3.10

# Validate .so
if [ ! -f "$SO_NAME" ]; then
    echo "✗ Python module build failed!"
    exit 1
fi

echo "✓ Python extension built successfully"


###############################################
# Final Summary
###############################################

echo ""
echo "=========================================="
echo "✓ Build Complete!"
echo "=========================================="
echo "Binary:        $CPP_REPO/bindings/openCOSMORS"
echo "Python module: $CPP_REPO/bindings/$SO_NAME"
echo ""
echo "You can test the Python module with:"
echo "  python3 run_example.py"
echo ""
echo "To rebuild:"
echo "  ./setup_scripts/build_opencosmo.sh $CPP_REPO"
echo "=========================================="
