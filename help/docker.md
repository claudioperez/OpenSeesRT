# Packaging and Cross Compiling

1. Build inside docker container of target OS
2. Build Python package (`python -m build`)
3. use [`auditwheel`](https://pypi.org/project/auditwheel/) to handle dynamic libs

```shell
{
docker run --rm dockcross/windows-shared-x86 
} > ./dockcross-windows-shared-x86
chmod +x ./dockcross-windows-shared-x86
```

> You may need to run `sudo dockerd` first to start
> the docker daemon.

Build numeric libraries:

```shell
cd OTHER && sudo ../dockcross_manylinux2014-x64 cmake -B bin_cross
sudo ../dockcross_manylinux2014-x64 bash -c 'cd bin_cross && make'
```

```shell
sudo ./cross bash # start bash session in container
sudo yum install tcl-devel
python3 -m pip install -U pip
python3 -m pip install pybind11
cmake .. -DDependencies=Cross -Dpybind11_ROOT=$(python3 -m pybind11 --cmakedir)
```


```
mkdir build_cross
cd build_cross
cp ../conanfile.py .
sudo ../dockcross-windows-shared-x86 conan install .
```

```shell
cd ../
sudo ./dockcross-windows-shared-x86 cmake -B build_win/
```


```
