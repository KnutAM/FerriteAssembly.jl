# External loading
In many simulations, there is a need to add contributions to the weak form 
in the form of external loads, either as source terms/body loads, or as 
non-zero Neumann boundary conditions. `FerriteAssembly` provides a 
`LoadHandler` for this purpose, such that it is not necessary to 
implement those terms into the specific elements directly. This avoid having 
to pass information to the elements about changes to the load levels etc. 

```@docs 
LoadHandler
Neumann
BodyLoad
DofLoad
``` 