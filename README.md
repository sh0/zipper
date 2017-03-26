## Mesh zippering

This repository contains 'mesh zippering' code, based on [ZipPack Polygon Mesh Zippering Package](https://graphics.stanford.edu/software/zippack/) version 1.14. The source has been ported from SGI IRIX system to be used on any modern operating system. User interface and mesh alignment code has been removed.

### Compiling and usage

There are no library dependencies, all you need is CMake build system and a C99 compatible compiler. To compile run:
```
mkdir build
cd build
cmake ..
make
```

Currently the program only zippers together two .ply format meshes using default parameters. The meshes should already contain triangle faces as no triangulation is carried out. Specify the mesh file paths as follows: `zipper src1.ply src2.ply dst.ply`.

Most of the original ZipPack functionality still exists, but is not enabled. If you need to change zippering parameters then inspect the constants in `zipper.c` file.

### Known issues

Some .ply files name vertex elements as `vertex_index`, but the zippering program expects `vertex_indices` identificator. If you get an error about reading .ply files, either change your mesh file header or change `#if 1` to `#if 0` on line 130 of `ply_wrapper.c`.

### License and references

The code is licensed under 3-clause BSD licence. Copyright belongs to Stanford University.

If you use this code in your research, please cite:
> "[Zippered Polygon Meshes from Range Images](https://graphics.stanford.edu/papers/zipper/zipper.pdf)"<br />
>  Greg Turk and Marc Levoy<br />
>  Proceedings of SIGGRAPH '94<br />
>  Pages 311-318
