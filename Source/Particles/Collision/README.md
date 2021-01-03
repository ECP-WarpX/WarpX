How to add a new collision operator:

Create a class that inherits from CollisionBase and include the method doCollisions
that does the operation. The constructor should call the CollisionBase constructor with
the collision name argument and read in any input parameters with a prefix of the
collision name that are specific to the operator.

Then modify CollisionHandler.cpp to include the header file of the new class.
In its constructor, add an if block for the input collision type creating an
instance of the new class in the allcollisions vector. Follow PairWiseCoulombCollision
as an example.
