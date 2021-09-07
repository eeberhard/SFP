# Spherical Frame Projection

Spherical frame projection (SFP) was conceived as a convenient
visualisation for orientational range of motion (ROM) data. The ROM
of a three-dimensional rotational joint, such as a ball joint, is
a subset of the group SO(3). If ROM is defined as the boundary between
feasible and infeasible orientation (pose), then it is the surface of a 
four-dimensional volume. This kind of thing is rather difficult to
represent visually or understand intuitively.

In most cases, Euler angles are used to reduce the dimensionality of the
orientation problem to three dimensions. However, a bounded ROM volume
generated from a sample of Euler angles may not accurately represent
or intuitively inform the reachability of joint due to inherent distortions
in the mapping, particularly near the poles.

The SFP approach projects the _frame_ of a pose onto the surface of a
unit sphere. The frame is represented by a set of three orthogonal basis
vectors, commonly labelled X, Y and Z, expressed in some base coordinate
system. If a pose is represented by a 3x3 homogeneous transformation
matrix, then the three projected frame points of that pose are just the
three columns of the matrix.

In a given sampling of ROM, a set of poses can be projected.
The result is a distribution of points on the surface of the sphere
which are divided into a triplet of regions; one for each frame axis.
The reachable area of each frame is represented by the respective points,
and a 2D surface boundary can be drawn around these points.

Each frame axis can then be imagined to be constrained by its respective
spherical surface boundary. Combining the constraints of each 2D surface
is enough to define and constrain the full rotational range.

This repository contains MATLAB code to generate and visualize an SFP
mapping from a given dataset of orientation samples.


TODO:
- Generate example SFP figure for README
- Add example folder with some ROM datasets
- Add simple usage instructions for SFP and SFPGUI

