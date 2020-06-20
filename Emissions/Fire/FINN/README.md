
# FINN Emissions

The executable, grid\_finn\_fire\_emis, is a fortran program that grids FINNv1.5 or FINNv2 fire emissions to the required input formats and grids for WRF-Chem, global regular lat-lon grids, and Spectral Element unstructured grids.  The same source code is used to produce any of the model formats, as specified in an input namelist.  Please address any questions about this code to the [WRF-Chem Fire\_Emis Forum](https://groups.google.com/a/ucar.edu/forum/?hl=en#!forum/wrf-chem-fire\_emiss)

Emissions files are available from [Christine Wiedinmyer](http://bai.acom.ucar.edu/Data/fire/)

*NOTE: FINNv1.5 and FINN2 files contain different species and different numbers of columns.  Please read the FINN2 README for more details (at: http://bai.acom.ucar.edu/Data/fire/).*

# Build
`cd src`
`make\_fire\_emis`

# Run
`cd run `
`../src/fire\_emis < {\_\_.inp}`

