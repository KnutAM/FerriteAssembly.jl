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
For brevity, the following abbrevations in type of analysis are used: 
DH=`DofHandler`, MDH=`MixedDofHandler` with `N` fieldhandlers, and 
MM = "Multiple materials on the grid, not connected to each fieldhandler"

For the levels, *Each X* refers to the datastructure of `CellBuffer` for the
cellset belonging to *X*, where *X* can be `MAT` (type of material), 
thread, or `FH` (fieldhandler `::FieldHandler` in `MixedDofHandler`).
Note that on the lower level (closer to the cell), the cellset is always 
reduced as the intersection between the higher level and the current. 
I.e. for multithreaded cases, *Each `FH`* refers to the cellset that 
is the intersection of the thread's cellset and the `FieldHandler`'s 
cellset. 


| Type of analysis  | Top level      | Each `MAT`   | Each thread  | Each `FH`    | Each cell    | 
| :---------------- | -------------- | ------------ | ------------ | ------------ | ------------ |
| DH                | `CellBuffer`   | -            | -            | -            | `CellBuffer` |
| MDH               | `Tuple`        | -            | -            | `CellBuffer` | `CellBuffer` |
| Threaded, DH      | `Vector`       | -            | `CellBuffer` | -            | `CellBuffer` |
| Threaded, MDH     | `Vector`       | -            | `Tuple`      | `CellBuffer` | `CellBuffer` |
| MM, DH            | `Dict{String}` | `CellBuffer` | -            | -            | `CellBuffer` |
| MM, MDH           | `Dict{String}` | `Tuple`      | -            | `CellBuffer` | `CellBuffer` |
| MM, Threaded, DH  | `Dict{String}` | `Vector`     | `CellBuffer` | -            | `CellBuffer` |
| MM, Threaded, MDH | `Dict{String}` | `Vector`     | `Tuple`      | `CellBuffer` | `CellBuffer` |

Note that the lowest level above *Each cell*, is always also the `CellBuffer`,
because one `CellBuffer` is shared between the cells that are being looped 
over in the innermost loop. 

## State variables 
`SV` is the datatype for one state variable. Normally, 
this exists for each integration point, and the lowest level (what is passed
into the `element_routine!`) is therefore a `Vector{SV}`, one `SV` per integration point.
It's datastructure is unaffected by threading: 
Each state entry belong to one cell so there are no race conditions. 

| Type of analysis  | Top level      | Each `MAT`   | Each `FH`    | Each cell    |
| :---------------- | -------------- | ------------ | ------------ | ------------ | 
| DH                | `Vector`       | -            | -            | `Vector{SV}` |
| MDH               | `Tuple`        | -            | `Dict{Int}`  | `Vector{SV}` |
| MM, DH            | `Dict{String}` | `Dict{Int}`  |              | `Vector{SV}` |
| MM, MDH           | `Dict{String}` | `Tuple`      | `Dict{Int}`  | `Vector{SV}` |

On the level above *Each cell*, whether that is a `Vector` or `Dict{Int}`, 
the indexing always refers to the global `cellid` (which is why a `Dict{Int}` is 
used when not all cells are included)

In the case of one `SV` per cell (when `create_states` are called without a cellvalue), 
the *Each cell* level becomes `SV` instead of `Vector{SV}` and the rest remain unchanged. 