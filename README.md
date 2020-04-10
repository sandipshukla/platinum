# platinum
Modern c++ path traced renderer.

This renderer is inspired from smallpt 
and uses obj scene files from Morgan Mcguire Free Models resoures.

![salle_de_bain_2048spp](https://github.com/sandipshukla/platinum/blob/master/salle_de_bain_2048spp.png)

Compile.

This has been developed on macOS catalina . To compile use the below command.

g++-9 -std=c++17 -o platinum -lstdc++fs -fopenmp platinum.cpp

Testing.

To test the renderer, first download the scenes folder from the release section.

https://github.com/sandipshukla/platinum/releases/download/v1.0/scenes.zip

unzip the scenes folder to the same location as the renderer. Run using the example command as given below.

platinum.out ./scenes/salle_de_bain/render.args

The renderer uses the render.args file to get settings for camera, image output, samples_per_pixel etc.

Each scene has its own render.args file.

