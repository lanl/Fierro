#!/bin/bash -e

# This script must be run from the top level Fierro directory
# This script assumes uncrustify has been installed



uncrustifyConfigDir="docs/formatting"
uncrustify_ex_dir="dev-utils/uncrustify/build"
# Save the directory the user had before calling
startDir=$(pwd)
thisScriptDir=$(dirname "$0")
echo "Start Directory: $startDir"
echo "This Script Directory: $thisScriptDir"
echo "Uncrustify config directory: $uncrustifyConfigDir"

# Get the full path to the config needed by uncrustify
cd $uncrustifyConfigDir
fullUncrustifyDir=$(pwd)
cd $startDir

#Get the full path to the uncrustify executable
cd $uncrustify_ex_dir
uncrustify_exe=$(pwd)/uncrustify
cd $startDir

# Get the source directory
cd src
sourceDir=$(pwd)

# Function to walk through files
treeProcess() {

  echo "calling tree in $(pwd)"
  # For every file in the directory
  for file in *.cpp; do
    if [ -f "$file" ]; then
      # echo "Using: $(pwd)/$f"
      "$uncrustify_exe" -c "$fullUncrustifyDir"/uncrustify.cfg --no-backup "$(pwd)/$file"
    fi
  done

  for file in *.h; do
    if [ -f "$file" ]; then
      # echo "Using: $(pwd)/$f"
      "$uncrustify_exe" -c "$fullUncrustifyDir"/uncrustify.cfg --no-backup "$(pwd)/$file"
    fi
  done

  # For every directory, recurse
  for dir in */; do
    if [ -d "$dir" ]; then
      cd "$dir" || exit
      treeProcess
      cd ..
    fi
    
  done
}


# Uncrustify SGH Solver
echo "Uncrusting: $sourceDir/Parallel-Solvers/Parallel-Explicit/SGH_Solver"
cd "$sourceDir/Parallel-Solvers/Parallel-Explicit/SGH_Solver" || exit
treeProcess

# Uncrustify Topology Optimization
echo "Uncrusting: $sourceDir/Parallel-Solvers/Parallel-Explicit/Topology_Optimization"
cd "$sourceDir/Parallel-Solvers/Parallel-Explicit/Topology_Optimization" || exit
treeProcess

# Uncrustify Eulerian Solver
echo "Uncrusting: $sourceDir/Parallel-Solvers/Parallel-Explicit/Eulerian_Solver"
cd "$sourceDir/Parallel-Solvers/Parallel-Explicit/Eulerian_Solver" || exit
treeProcess

# Uncrustify Dynamic Elastic
echo "Uncrusting: $sourceDir/Parallel-Solvers/Parallel-Explicit/Dynamic_Elastic_Solver"
cd "$sourceDir/Parallel-Solvers/Parallel-Explicit/Dynamic_Elastic_Solver" || exit
treeProcess

exit

