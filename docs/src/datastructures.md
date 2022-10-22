# Data structures
The assembly requires different datastructure, which also differ depending on 
how the analysis is setup. An important feature of this package is that it provides
a suggestion for how these structures can be organized to work in a range of cases. 
While it is possible to use the package without knowing exactly how these are made, 
it will make tracking down bugs easier and may be important for postprocessing.

## CellBuffer
All values passed into `element_routine!`, except the state variables, 
`dh_fh::Union{DofHandler,FieldHandler}`, and time increment,
are passed taken from the relevant `CellBuffer`. 
(The cellbuffer itself is also passed into `element_routine!`. 
By looking at the created `CellBuffer`, it is therefore clear what will 
be passed into your `element_routine!`. 

The following items deserve extra attention

* `cellvalues`: Can be either one `CellValues` (e.g `CellVectorValues`), 
  or a collection of these. This collection must be either a `Tuple` or 
  a `NamedTuple` for the reinitialization will work. If provided as a custom 
  type, just define `Ferrite.reinit!(cv::MyCustomCellValueCollection, coords)`
* `material`: Represents the material to dispatch on for that specific element. 
* `cell_load`: Is an additional custom user type, intended for defining source 
  terms or body loads. 
* `cache`: A final custom user type, intended to store preallocated data 
  (e.g. for `FE²` simulations)


Depending on the analysis, multiple `CellBuffers` may be created and 
nested in various structures. These structures are explained in 
the following table. 
For brevity, the following abbrevations are used: 
DH=`DofHandler`, MDH=`MixedDofHandler` with `N` fieldhandlers
MM = "Multiple materials on the grid, not connected to each fieldhandler"

| Type of analysis  | Top level (TL) | TL-1         | TL-2         | TL-3         |
| :---------------- | -------------- | ------------ | ------------ | -----------  | 
| DH                | `CellBuffer`   |              |              |              |
| MDH               | `Tuple`        | `CellBuffer` |              |              |    
| Threaded, DH      | `Vector`       | `CellBuffer` |              |              | 
| Threaded, MDH     | `Vector`       | `Tuple`      | `CellBuffer` |              | 
| DH, MM            | `Dict{String}` | `CellBuffer` |              |              | 
| MDH, MM           | `Dict{String}` | `Tuple`      | `CellBuffer` |              | 
| Threaded, DH, MM  | `Dict{String}` | `Vector`     | `CellBuffer` |              | 
| Threaded, MDH, MM | `Dict{String}` | `Vector`     | `Tuple`      | `CellBuffer` |


## State variables 
`SV` is the datatype for one state variable. Normally, 
this exists for each integration point, and the lowest level (what is passed
into the `element_routine!`) is therefore a `Vector{SV}`, one `SV` per integration point.
It's datastructure is unaffected by threading: 
Each state entry belong to one cell so there are no race conditions. 

| Type of analysis  | Top level (TL) | TL-1         | TL-2         | TL-3         |
| :---------------- | -------------- | ------------ | ------------ | ------------ | 
| DH                | `Vector`       | `Vector{SV}` |              |              |
| MDH               | `Tuple`        | `Dict{Int}`  | `Vector{SV}` |              |
| DH, MM            | `Dict{String}` | `Dict{Int}`  | `Vector{SV}` |              |
| MDH, MM           | `Dict{String}` | `Tuple`      | `Dict{Int}`  | `Vector{SV}` |

In the case of one `SV` per cell when `create_states` are called without a cellvalue, 
the lowest level becomes `SV` instead of `Vector{SV}`
