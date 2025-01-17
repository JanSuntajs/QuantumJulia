using DrWatson
@quickactivate :QuantumJulia

using Literate
using Printf

repo_src = projectdir()
notebooks_dir = projectdir("notebooks")

function build_notebooks_from_jl(folder::String=notebooks_dir; output_folder::String = projectdir("generated_notebooks"),
    execute::Bool=false)

    succes_count = 0
    failure_count = 0

    # Create the output folder if it doesn't exist
    isdir(output_folder) || mkdir(output_folder)

    # Get a list of all .jl files in the specified folder
    jl_files = filter(f -> endswith(f, ".jl"), readdir(folder, join=true))

    for jl_file in jl_files
        # Define the output .ipynb file name
        notebook_name = joinpath(output_folder, basename(replace(jl_file, ".jl" => ".ipynb")))
        
        # Convert the .jl file to a .ipynb notebook
        println("Building $jl_file into $notebook_name...")
        try
            Literate.notebook(jl_file, output_folder; execute=execute)
            println("Successfully built $notebook_name.")
            succes_count += 1
        catch e
            println("Failed to build $jl_file: $e")
            failure_count += 1
        end
    end

    if failure_count == 0
        println("All .jl files in $folder have been converted to .ipynb files in $output_folder.")
    else
        println("Successfully converted $success_count files.")
        println("Failed to convert $failure_count files.")
    end

end

build_notebooks_from_jl()