#!/bin/bash

# SOURCE: 
#  https://giswqs.medium.com/building-a-conda-package-and-uploading-it-to-anaconda-cloud-6a3abd1c5c52

# change the package name to the existing PyPi package you would like to build
# and adjust the Python versions

pkg='opensees'
recipe="etc/conda"
python_versions=( 3.7 3.8 3.9 )
platforms=( osx-64 linux-32 linux-64 win-32 win-64 )

alias CONDA-BUILD="conda mambabuild"

if "" then
  echo "Building conda package ..."
  cd ~
  conda skeleton pypi $pkg
  cd $pkg
  wget https://conda.io/docs/_downloads/build1.sh
  wget https://conda.io/docs/_downloads/bld.bat
  cd ~
fi

# building conda packages
for v in "${python_versions[@]}"
do
	CONDA-BUILD --python $v $recipe
done

# convert package to other platforms
cd ~
find $HOME/*/conda-bld/linux-64/ -name *.tar.bz2 | while read file
do
    echo "converting $file"
    #conda convert --platform all $file  -o $HOME/conda-bld/
    for platform in "${platforms[@]}"
    do
       conda convert --platform $platform $file  -o $HOME/conda-bld/
    done
    
done

# upload packages to conda
find $HOME/conda-bld/ -name *.tar.bz2 | while read file
do
    echo "uploading $file"
    anaconda upload $file
done

echo "Building conda package complete!"
