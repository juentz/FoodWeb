"
Make sure Excel files are ordered in the right way!
Takes information out of the provided excel files in the matrices folder
Defines a list of all species with its properties
Defines different communties that shall be analyzed with its properties

Possibility to put some extra knowledge in a single Community

saves communities in extra file -name- which is called by -filenameApplication-
"

using XLSX # for extracting the matrices out of excel files

"GET DATA"
# Open Excel file
xf = XLSX.readxlsx("matrices/Arctic.xlsx")
# Get worksheet
sheet = xf["PredPreyMatrixArctic1"]
# Extract foodweb
# A = Float64.(coalesce.(XLSX.getdata(sheet, "D4:BS71"), 0.0))
