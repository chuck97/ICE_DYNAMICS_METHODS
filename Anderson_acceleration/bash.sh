cd src
rm *.o *.exe *.mod
echo 1
cd ..
cd bin
rm *.o *.exe *.mod *.ps
echo 2
cd ..
cd obj
rm *.o *.exe *.mod
echo 3
cd ..
cd mod
rm *.o *.exe *.mod
echo 4
cd ..
cd graphics 
rm *.vtk
echo 5
cd ..
cd resudal
rm *.txt
echo 6
cd ..
cd values
cd nodals 
rm *.txt
echo 7
cd ..
cd triangles
rm *.txt
echo 8
cd ..
cd ..
cd src
gfortran -fPIC -g -Wline-truncation -fcheck=all -c Classes.f95 Constants.f95 Numerical_integration.f95 Grid.f95 Assembling.f95 Matricies.f95 Flux_correction.f95 Dynamics.f95 Main.f95
echo 9
mv *.o ../obj/
mv *.mod ../mod/
rm *mod0
cd ..
cd obj
gfortran -fPIC Classes.o Constants.o Numerical_integration.o Grid.o Assembling.o Matricies.o Flux_correction.o Dynamics.o Main.o -L ../../Libraries/ -lblas-3.2 -llapack-3.2 -llapack_ext-3.2 -llu-5.1  -lilu-3.1 -lblas-3.2 -llapack-3.2 -llapack_ext-3.2 -laft2D-3.1 -lview2D-3.1 -lmba2D-3.1 -lskit -llapack -lrefblas -ltmglib
echo 10
mv a.out result.exe
mv result.exe ../bin/
cd ..
cd bin
./result.exe
 
