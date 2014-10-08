#MeshValmet: Validation Metric for Meshes

##What is it?

MeshValmet is a tool that measures surface to surface distance between two triangle meshes using user-specified uniform sampling. Thus, users can choose finer sampling level to calculate errors to gain more accuracy in the "error space", or sparser sampling to gain speed and get an approximate feeling of error distribution between boundaries.

Besides its pleasant visualization using the VTK library, MeshValmet also provides useful histogram and statistical information based on the sample errors, such as mean and median distance, root mean square distance, mean square distance, mean absolute distance, Hausdorff distance, 95 percentile, 68 percentile, etc.

MeshValmet is based on the work of Nicolas Aspert, etc.: MESH: Measuring Errors between Surfaces using the Hausdorff distance in the proceedings of the IEEE Int. Conf. on Multimedia and Expo 2002 (ICME), vol. I, pp. 705-708.

The calculation of the Dice's Coefficient is calculated by Joshua Stough using the concept of a Riemannian sum.

##License

See License.txt

##More information

Find the tool on [NITRC](http://www.nitrc.org/projects/meshvalmet)

