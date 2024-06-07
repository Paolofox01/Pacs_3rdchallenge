### Third challenge of the APSC course 

## Content

`laplace_solver.hpp` -> definition of a class that solves the Laplace equation

`laplace_solver.cpp` -> implementaation of the members of laplace_solver.hpp

`main.cpp` -> initializes the grid and launches the solver

`parameters.json` -> contains some parameters to use for the algorithm

`Doxifile` -> contains the instructions needed to create the documentation using Doxygen

`Makefile` -> contains the intructions needed to compile the code

## How to compile and run the code

modify the PACS_ROOT path in the Makefile and then run 

````
make
````

to run the code using the matrix written in the json file just run

````
mpiexec -n i ./main j k
````
where i is the number of processors to use with MPI, j the length of the side of the grid and k the openMP cores,
j and k can also be changed from the json file



if you want you can also run 

````
./test/test_scalability.sh
````
for some examples with chosen numbers for the parameters, the results will be stored in the data directory


the code creates a .vtk with the solution. You can visualizze it with 

````
paraview filename.vtk
````

## How to create the documentiation with Doxygen

run 

````
make doc
````

then open the index.html file in /doc/html to view the documentation 

## disclaimer

I tried passing the function with muparser but I had some troubles, I left my attempt in the code commmented and I just used a function for f(x)

