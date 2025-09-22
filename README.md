# Mandel.jl

This is a Julia package integrated with [Pluto.jl](https://plutojl.org/) to visualize complex dynamical systems.

At this moment, there are two user interfaces **MandelPluto** and **MandelMakie**. The first one is embedded directly on Pluto notebooks, but it is still in the early stages of development, and the second one uses the **GLMakie** backend of [**Makie.jl**](https://docs.makie.org/stable/) and opens a separate window to interact with plots.

## How to use MandelMakie

### If you already have Pluto.jl installed

First, run Julia with multithreading enabled
```
julia -t auto
```
Then, run Pluto
```julia
import Pluto
Pluto.run()
```
Lastly, run the following code on a new **Pluto.jl** notebook.
```julia
begin
    import Pkg
    Pkg.activate(mktempdir())
    Pkg.add(url="https://github.com/vambrosi/Mandel.jl.git",
			subdir="src/MandelMakie")
    using MandelMakie
end
```
This will activate a temporary environment, but you can also [use your on environment](https://plutojl.org/en/docs/packages-advanced/).

You can then write `?Viewer` in a cell to see some examples and learn about the user interface.

### Using the REPL

You can also run Julia as above, but add **MandelMakie** directly in the REPL:
```julia
import Pkg
Pkg.activate(mktempdir())
Pkg.add(url="https://github.com/vambrosi/Mandel.jl.git",
  subdir="src/MandelMakie")
using MandelMakie
```
Running `?Viewer` in the REPL will show you the documentation.

### Using Git

Clone this repository, and go to the folder `./src/MandelMakie` on your terminal.

Run Julia as before, then activate and instantiate the **MandelMakie** environment
```julia
]activate .
instantiate
```
You will only need to run `instantiate` once to install all the dependencies, but you will need to run `activate` every time you want to use this package. The `.` after `activate` represents the current folder where the `Project.toml` file with the dependencies is located.

Press `Ctrl+C` or `backspace` to exit the package manager and run the following lines to start Pluto
```julia
import Pluto
Pluto.run()
```
This will open a browser window, and you can open the file `Examples.jl` for further instructions and examples of how to interact with the package.

### Using Visual Studio Code

This is similar to the last section, so we will mostly point out the differences.

To use multithreading, click on *Extensions -> Julia -> &#9881; -> Extension Settings* and edit **Julia: Num Threads** to `"auto"`.

Open the folder `./src/MandelMakie` on [**Visual Studio Code**](https://code.visualstudio.com/docs/languages/julia).

Press `Alt + J, Alt + O` to open a Julia REPL. By default VSCode will `activate` the environment in the folder (you will see `Julia env: MandelMakie` in the bottom of the window). So just run `instantiate` the first time you use the package. Afterwards, you will only need to `import` and `run` Pluto.
